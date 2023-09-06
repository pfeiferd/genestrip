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
import java.io.PrintStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.GSConfig;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.match.FastqKMerMatcher;
import org.metagene.genestrip.match.FastqKMerMatcher.Result;
import org.metagene.genestrip.match.ResultReporter;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.store.KMerUniqueCounter;
import org.metagene.genestrip.store.KMerUniqueCounterBits;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;

public class MultiMatchGoal extends FileListGoal<GSProject> {
	public static final String NAME = "multimatch";
	
	private final Map<File, List<File>> fileToFastqs;
	private final File csvFile;
	private final KMerStoreFileGoal storeGoal;
	private final boolean writedFiltered;

	private FastqKMerMatcher matcher;
	private KMerStoreWrapper wrapper;
	private ResultReporter reporter;
	private KMerUniqueCounter uniqueCounter;

	@SafeVarargs
	public MultiMatchGoal(GSProject project, String name, File csvFile, KMerStoreFileGoal storeGoal,
			boolean writeFiltered, Goal<GSProject>... deps) {
		super(project, name, (List<File>) null, ArraysUtil.append(deps, storeGoal));
		this.csvFile = csvFile;
		this.storeGoal = storeGoal;
		this.writedFiltered = writeFiltered;
		fileToFastqs = new HashMap<File, List<File>>();
	}

	@Override
	protected void provideFiles() {
		try {
			CSVParser parser = readCSVFile(csvFile);

			for (CSVRecord record : parser) {
				String prefix = record.get(0);
				String fastqFilePath = record.get(1);
				File fastq = new File(fastqFilePath);
				if (fastq.exists()) {
					File matchFile = getProject().getOutputFile(prefix + "_" + getName(), null, FileType.CSV, false);
					List<File> fastqs = fileToFastqs.get(matchFile);
					if (fastqs == null) {
						fastqs = new ArrayList<File>();
						fileToFastqs.put(matchFile, fastqs);
						addFile(matchFile);
					}
					fastqs.add(fastq);
				} else {
					getLogger().warn("Ignoring missing fastq file " + fastq + ".");
				}
			}

			parser.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	
	@Override
	protected void makeFile(File file) {
		try {
			File filteredFile = null;
			File krakenOutStyleFile = null;
			List<File> fastqs = fileToFastqs.get(file);

			if (writedFiltered) {
				filteredFile = getProject().getOutputFile(getName(), file, FileType.FASTQ_RES, true);
				krakenOutStyleFile = getProject().getOutputFile(getName(), file, FileType.KRAKEN_OUT_RES, false);
			}

			GSConfig config = getProject().getConfig();
			if (matcher == null) {
				wrapper = KMerStoreWrapper.load(storeGoal.getFile());
				matcher = new FastqKMerMatcher(wrapper.getKmerStore(), config.getMaxReadSizeBytes(),
						config.getThreadQueueSize(), config.getThreads());
				reporter = new ResultReporter(wrapper.getTaxids());
				uniqueCounter = config.isCountUniqueKmers() ? new KMerUniqueCounterBits(wrapper.getKmerStore()) : null;
			}
			if (uniqueCounter != null) {
				uniqueCounter.clear();
			}
			Result res = matcher.runClassifier(fastqs, filteredFile, krakenOutStyleFile, uniqueCounter);
			PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
			reporter.printMatchResult(res, out, wrapper.getStoreStats());
			out.close();
		} catch (IOException | ClassNotFoundException e) {
			throw new RuntimeException(e);
		}
	}
	
	public static CSVParser readCSVFile(File csvFile) throws IOException {
		Reader in = new InputStreamReader(StreamProvider.getInputStreamForFile(csvFile));
		CSVFormat format = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#').setDelimiter(' ')
				.setRecordSeparator('\n').build();
		return format.parse(in);		
	}

	@Override
	protected void endMake() {
		if (matcher != null) {
			matcher.dump();
			matcher = null;
			wrapper = null;
			uniqueCounter = null;
			reporter = null;
		}
	}

	protected static class PrefixAndFile {
		protected String prefix;
		protected File file;

		public PrefixAndFile(String prefix, File file) {
			this.prefix = prefix;
			this.file = file;
		}
	}
}
