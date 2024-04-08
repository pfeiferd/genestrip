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
package org.metagene.genestrip.goals;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.logging.Log;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;

public abstract class MultiFileGoal extends FileListGoal<GSProject> {
	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#')
			.setDelimiter(' ').setRecordSeparator('\n').build();

	private final boolean csv;
	private final File csvOrFastqFile;
	protected final Map<File, List<StreamingResource>> fileToFastqs;

	@SafeVarargs
	public MultiFileGoal(GSProject project, String name, Goal<GSProject>... deps) {
		super(project, name, (List<File>) null, deps);
		this.csv = true;
		this.csvOrFastqFile = null;
		fileToFastqs = new HashMap<File, List<StreamingResource>>();
	}

	@SafeVarargs
	public MultiFileGoal(GSProject project, String name, boolean csv, File csvOrFastqFile, Goal<GSProject>... deps) {
		super(project, name, (List<File>) null, deps);
		this.csv = csv;
		this.csvOrFastqFile = csvOrFastqFile;
		fileToFastqs = new HashMap<File, List<StreamingResource>>();
	}

	protected File getSourceDir() {
		return getProject().getFastqDir();
	}

	@Override
	protected void provideFiles() {
		Map<String, List<StreamingResource>> keyToFastqs = getProject().getKeyToFastqs();
		if (keyToFastqs == null) {
			if (csv) {
				keyToFastqs = readMultiCSV(getSourceDir(), csvOrFastqFile, getLogger());
			} else {
				keyToFastqs = new HashMap<String, List<StreamingResource>>();
				keyToFastqs.put(null, Collections.singletonList(new StreamingFileResource(csvOrFastqFile)));
			}
		}
		for (String key : keyToFastqs.keySet()) {
			File matchFile;
			if (key != null) {
				matchFile = getProject().getOutputFile(getName(), key, null, FileType.CSV, false);
			} else {
				StreamingResource file = keyToFastqs.get(null).get(0);
				matchFile = getProject().getOutputFile(getName(), null, file.getName(), FileType.CSV, false);
			}
			addFile(matchFile);
			fileToFastqs.put(matchFile, keyToFastqs.get(key));
		}
	}

	// We use a linked hash map because it preserves the order of the keys from the
	// file.
	public static LinkedHashMap<String, List<StreamingResource>> readMultiCSV(File defaultDir, File csvFile, Log logger) {
		try {
			LinkedHashMap<String, List<StreamingResource>> res = new LinkedHashMap<String, List<StreamingResource>>();
			CSVParser parser;
			parser = readCSVFile(csvFile);
			for (CSVRecord record : parser) {
				String name = record.get(0);
				String fastqFilePath = record.get(1);
				File fastq = new File(fastqFilePath);
				if (!fastq.exists()) {
					fastq = new File(defaultDir, fastqFilePath);
				}
				if (fastq.exists()) {
					List<StreamingResource> fastqs = res.get(name);
					if (fastqs == null) {
						fastqs = new ArrayList<StreamingResource>();
						res.put(name, fastqs);
					}
					fastqs.add(new StreamingFileResource(fastq));
				} else {
					if (logger.isWarnEnabled()) {
						logger.warn("Ignoring missing file " + fastq + ".");
					}
				}
			}

			parser.close();

			return res;
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public static CSVParser readCSVFile(File csvFile) throws IOException {
		Reader in = new InputStreamReader(StreamProvider.getInputStreamForFile(csvFile));
		return FORMAT.parse(in);
	}
}
