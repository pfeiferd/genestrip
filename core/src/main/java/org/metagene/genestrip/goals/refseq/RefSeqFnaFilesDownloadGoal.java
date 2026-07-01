/*
 * 
 * “Commons Clause” License Condition v1.0
 * 
 * The Software is provided to you by the Licensor under the License, 
 * as defined below, subject to the following condition.
 * 
 * Without limiting other conditions in the License, the grant of rights under the License 
 * will not include, and the License does not grant to you, the right to Sell the Software.
 * 
 * For purposes of the foregoing, “Sell” means practicing any or all of the rights granted 
 * to you under the License to provide to third parties, for a fee or other consideration 
 * (including without limitation fees for hosting or consulting/ support services related to 
 * the Software), a product or service whose value derives, entirely or substantially, from the 
 * functionality of the Software. Any license notice or attribution required by the License 
 * must also include this Commons Clause License Condition notice.
 * 
 * Software: genestrip
 * 
 * License: Apache 2.0
 * 
 * Licensor: Daniel Pfeifer (daniel.pfeifer@progotec.de)
 * 
 */
package org.metagene.genestrip.goals.refseq;
import java.nio.charset.StandardCharsets;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSConfigKey.SeqType;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.RefSeqCategory;

/**
 * Download goal for the RefSeq genomic/RNA {@code .fna} files of the selected categories, verified
 * against their catalog MD5 checksums; the list of files to fetch is derived from the release
 * {@code files.installed} listing and the configured sequence type.
 *
 * @param <P> the project type
 */
public class RefSeqFnaFilesDownloadGoal<P extends GSProject> extends RefSeqDownloadGoal<P> {
	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setDelimiter('\t').setRecordSeparator('\n')
			.build();

	private final ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal;
	private final RefSeqCatalogDownloadGoal catalogDLGoal;
	private final ObjectGoal<CheckRefSeqRNumGoal.Result, P> checkReleaseNGoal;
	private List<File> files;
	private Map<File, RefSeqCategory> file2Cat;
	private ObjectGoal<Map<String, String>, P> checkSumGoal;

	/**
	 * Creates the goal, wiring the categories, catalog download, checksum and release-check goals.
	 *
	 * @param project           the project this goal belongs to
	 * @param categoriesGoal    the goal supplying the RefSeq categories to download
	 * @param catalogDLGoal     the goal downloading the RefSeq catalog
	 * @param checkSumGoal      the goal supplying the expected file checksums
	 * @param checkReleaseNGoal the goal checking the RefSeq release number
	 * @param deps              any further goals this goal depends on
	 */
	@SafeVarargs
	public RefSeqFnaFilesDownloadGoal(P project,
			ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal, RefSeqCatalogDownloadGoal catalogDLGoal,
			ObjectGoal<Map<String, String>, P> checkSumGoal, ObjectGoal<CheckRefSeqRNumGoal.Result, P> checkReleaseNGoal,
			Goal<P>... deps) {
		super(project, GSGoalKey.REFSEQFNA, Goal.append(deps, categoriesGoal, catalogDLGoal, checkSumGoal, checkReleaseNGoal));

		this.categoriesGoal = categoriesGoal;
		this.catalogDLGoal = catalogDLGoal;
		this.checkSumGoal = checkSumGoal;
		this.checkReleaseNGoal = checkReleaseNGoal;
	}

	/**
	 * Runs the RefSeq release check before reporting whether this goal is already made.
	 */
	@Override
	public boolean isMade() {
		checkReleaseNGoal.get();
		return super.isMade();
	}

	/**
	 * Downloads the given RefSeq FASTA file, first warning if the local release is outdated (in which
	 * case the MD5 check may fail).
	 */
	@Override
	protected void makeFile(File file) throws IOException {
		if (CheckRefSeqRNumGoal.Result.OUTDATED.equals(checkReleaseNGoal.get())) {
			if (getLogger().isWarnEnabled()) {
				getLogger().warn("MD5 check will likely fail as you are working with an outdated RefSeq release.");
			}
		}

		super.makeFile(file);
	}

	/**
	 * The RefSeq category that the given downloaded file belongs to.
	 *
	 * @param file the downloaded file
	 * @return the category the file belongs to, or {@code null} if unknown
	 */
	public RefSeqCategory getCategoryForFile(File file) {
		return file2Cat.get(file);
	}

	/**
	 * Whether the given catalog file name matches the configured sequence type (genomic and/or RNA).
	 *
	 * @param filename the catalog file name to check
	 * @return {@code true} if the file name matches the configured sequence type
	 */
	protected boolean isRelevantFileName(String filename) {
		// TODO: Code not very elegant - but whatever...
		switch ((SeqType) configValue(GSConfigKey.SEQ_TYPE)) {
		case RNA:
		case M_RNA:
		case ALL_RNA:
			return filename.endsWith(".rna.fna.gz") || filename.endsWith(".rna.fna");
		case ALL:
			return filename.endsWith(".genomic.fna.gz") || filename.endsWith(".genomic.fna")
					|| filename.endsWith(".rna.fna.gz") || filename.endsWith(".rna.fna");
		case GENOMIC:
		default:
			return filename.endsWith(".genomic.fna.gz") || filename.endsWith(".genomic.fna");
		}
	}

	/**
	 * The selected {@link RefSeqCategory} whose directory prefixes the given file name, or
	 * {@code null} if none matches.
	 *
	 * @param filename the file name to match
	 * @return the matching category, or {@code null} if none matches
	 */
	protected RefSeqCategory getCategoryForFileName(String filename) {
		for (RefSeqCategory cat : categoriesGoal.get()) {
			if (filename.startsWith(cat.getDirectory() + ".")) {
				return cat;
			}
		}
		return null;
	}

	@Override
	protected String getFTPDir(File file) {
		return RELEASE_FOLDER + "/" + file2Cat.get(file).getDirectory();
	}
	
	@Override
	protected String getMD5CheckSum(File file) {
		return checkSumGoal.get().get(file.getName());
	}

	/**
	 * Builds the deterministically ordered list of {@code .fna} files to download by scanning the
	 * catalog's {@code files.installed} for relevant files in the selected categories.
	 */
	@Override
	public List<File> getFiles() {
		if (files == null) {
			// There is no other way: The catalog file must be there early, in order to get
			// the file list....
			catalogDLGoal.make();

			files = new ArrayList<File>();
			file2Cat = new HashMap<File, RefSeqCategory>();

			File installedFiles = catalogDLGoal.getInstalledFilesFile();

			try (Reader in = new InputStreamReader(StreamProvider.getInputStreamForFile(installedFiles), StandardCharsets.UTF_8)) {
				Iterable<CSVRecord> records = FORMAT.parse(in);
				for (CSVRecord record : records) {
					String filename = record.get(1);
					RefSeqCategory cat = getCategoryForFileName(filename);
					if (cat != null) {
						if (isRelevantFileName(filename)) {
							File file = new File(getProject().getCommon().getRefSeqDir(), filename);
							files.add(file);
							file2Cat.put(file, cat);
						}
					}
				}
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
			// It is very important to get the files in a well determined order for the
			// later steps.
			// (Any order would be fine.)
			Collections.sort(files);
		}
		return files;
	}
}
