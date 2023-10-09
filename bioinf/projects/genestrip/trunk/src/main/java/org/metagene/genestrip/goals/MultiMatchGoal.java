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
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.FastqKMerMatcher;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.match.ResultReporter;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.store.KMerUniqueCounter;
import org.metagene.genestrip.store.KMerUniqueCounterBits;
import org.metagene.genestrip.tax.TaxTree;

public class MultiMatchGoal extends FileListGoal<GSProject> {
	public static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#')
			.setDelimiter(' ').setRecordSeparator('\n').build();

	public static final String NAME = "multimatch";

	private final Map<File, List<File>> fileToFastqs;
	private final File csvFile;
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final ObjectGoal<KMerStoreWrapper, GSProject> storeGoal;
	private final boolean writedFiltered;

	private FastqKMerMatcher matcher;
	private KMerStoreWrapper wrapper;
	private ResultReporter reporter;
	private KMerUniqueCounter uniqueCounter;

	@SafeVarargs
	public MultiMatchGoal(GSProject project, String name, File csvFile, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			ObjectGoal<KMerStoreWrapper, GSProject> storeGoal, boolean writeFiltered, Goal<GSProject>... deps) {
		super(project, name, (List<File>) null, Goal.append(deps, taxTreeGoal, storeGoal));
		this.csvFile = csvFile;
		this.taxTreeGoal = taxTreeGoal;
		this.storeGoal = storeGoal;
		this.writedFiltered = writeFiltered;
		fileToFastqs = new HashMap<File, List<File>>();
	}

	@Override
	protected void provideFiles() {
		try {
			CSVParser parser = readCSVFile(csvFile);

			for (CSVRecord record : parser) {
				String name = record.get(0);
				String fastqFilePath = record.get(1);
				File fastq = new File(fastqFilePath);
				if (!fastq.exists()) {
					fastq = new File(getProject().getFastqDir(), fastqFilePath);
				}
				if (fastq.exists()) {
					File matchFile = getProject().getOutputFile(getName() + "_" + name, null, FileType.CSV, false);
					List<File> fastqs = fileToFastqs.get(matchFile);
					if (fastqs == null) {
						fastqs = new ArrayList<File>();
						fileToFastqs.put(matchFile, fastqs);
						addFile(matchFile);
					}
					fastqs.add(fastq);
				} else {
					if (getLogger().isWarnEnabled()) {
						getLogger().warn("Ignoring missing fastq file " + fastq + ".");
					}
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
				wrapper = storeGoal.get();
				matcher = new FastqKMerMatcher(wrapper.getKmerStore(), config.getMaxReadSizeBytes(),
						config.getThreadQueueSize(), config.getThreads(), config.getMaxKMerResCounts());
				reporter = new ResultReporter(taxTreeGoal.get());
				uniqueCounter = config.isCountUniqueKMers() ? new KMerUniqueCounterBits(wrapper.getKmerStore(), true)
						: null;
			}
			if (uniqueCounter != null) {
				uniqueCounter.clear();
			}
			MatchingResult res = matcher.runClassifier(fastqs, filteredFile, krakenOutStyleFile, uniqueCounter);
			PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
			reporter.printMatchResult(res, out, wrapper);
			out.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public static CSVParser readCSVFile(File csvFile) throws IOException {
		Reader in = new InputStreamReader(StreamProvider.getInputStreamForFile(csvFile));
		return FORMAT.parse(in);
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
