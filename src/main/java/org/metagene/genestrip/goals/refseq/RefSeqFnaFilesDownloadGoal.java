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
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;

public class RefSeqFnaFilesDownloadGoal extends RefSeqDownloadGoal {
	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setDelimiter('\t').setRecordSeparator('\n')
			.build();

	private final ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal;
	private final RefSeqCatalogDownloadGoal catalogDLGoal;
	private List<File> files;
	private Map<File, RefSeqCategory> file2Cat;

	@SafeVarargs
	public RefSeqFnaFilesDownloadGoal(GSProject project, String name,
			ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal, RefSeqCatalogDownloadGoal catalogDLGoal,
			Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, categoriesGoal, catalogDLGoal));

		this.categoriesGoal = categoriesGoal;
		this.catalogDLGoal = catalogDLGoal;
	}

	public RefSeqCategory getCategoryForFile(File file) {
		return file2Cat.get(file);
	}

	protected boolean isRelevantFileName(String filename) {
		return filename.endsWith(".genomic.fna.gz") || filename.endsWith(".genomic.fna");
	}

	protected RefSeqCategory getCategoryForFileName(String filename) {
		for (RefSeqCategory cat : categoriesGoal.get()[1]) {
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
	public List<File> getFiles() {
		if (files == null) {
			// There is no other way: The catalog file must be there early, in order to get
			// the file list....
			catalogDLGoal.make();

			files = new ArrayList<File>();
			file2Cat = new HashMap<File, RefSeqCategory>();

			File installedFiles = catalogDLGoal.getInstalledFilesFile();

			try {
				Reader in = new InputStreamReader(StreamProvider.getInputStreamForFile(installedFiles));
				Iterable<CSVRecord> records = FORMAT.parse(in);
				for (CSVRecord record : records) {
					String filename = record.get(1);
					RefSeqCategory cat = getCategoryForFileName(filename);
					if (cat != null) {
						if (isRelevantFileName(filename)) {
							File file = new File(getProject().getConfig().getRefSeqDir(), filename);
							files.add(file);
							file2Cat.put(file, cat);
						}
					}
				}
				in.close();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
			Collections.sort(files);
		}
		return files;
	}
}
