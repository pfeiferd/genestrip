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

public class RefSeqFnaFilesDownloadGoal extends RefSeqDownloadGoal {
	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setDelimiter('\t').setRecordSeparator('\n')
			.build();

	private final ObjectGoal<Set<RefSeqCategory>, GSProject> categoriesGoal;
	private final RefSeqCatalogDownloadGoal catalogDLGoal;
	private final ObjectGoal<CheckRefSeqRNumGoal.Result, GSProject> checkReleaseNGoal;
	private List<File> files;
	private Map<File, RefSeqCategory> file2Cat;
	private ObjectGoal<Map<String, String>, GSProject> checkSumGoal;

	@SafeVarargs
	public RefSeqFnaFilesDownloadGoal(GSProject project, 
			ObjectGoal<Set<RefSeqCategory>, GSProject> categoriesGoal, RefSeqCatalogDownloadGoal catalogDLGoal,
			ObjectGoal<Map<String, String>, GSProject> checkSumGoal, ObjectGoal<CheckRefSeqRNumGoal.Result, GSProject> checkReleaseNGoal,
			Goal<GSProject>... deps) {
		super(project, GSGoalKey.REFSEQFNA, Goal.append(deps, categoriesGoal, catalogDLGoal, checkSumGoal, checkReleaseNGoal));

		this.categoriesGoal = categoriesGoal;
		this.catalogDLGoal = catalogDLGoal;
		this.checkSumGoal = checkSumGoal;
		this.checkReleaseNGoal = checkReleaseNGoal;
	}

	@Override
	public boolean isMade() {
		checkReleaseNGoal.get();
		return super.isMade();
	}

	@Override
	protected void makeFile(File file) throws IOException {
		if (CheckRefSeqRNumGoal.Result.OUTDATED.equals(checkReleaseNGoal.get())) {
			if (getLogger().isWarnEnabled()) {
				getLogger().warn("MD5 check will likely fail as you are working with an outdated RefSeq release.");
			}
		}

		super.makeFile(file);
	}

	public RefSeqCategory getCategoryForFile(File file) {
		return file2Cat.get(file);
	}

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

	@Override
	public List<File> getFiles() {
		if (files == null) {
			// There is no other way: The catalog file must be there early, in order to get
			// the file list....
			catalogDLGoal.make();

			files = new ArrayList<File>();
			file2Cat = new HashMap<File, RefSeqCategory>();

			File installedFiles = catalogDLGoal.getInstalledFilesFile();

			try (Reader in = new InputStreamReader(StreamProvider.getInputStreamForFile(installedFiles))) {
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
