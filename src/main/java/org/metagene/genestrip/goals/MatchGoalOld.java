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
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.FastqKMerMatcher;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.match.ResultReporter;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerUniqueCounter;
import org.metagene.genestrip.store.KMerUniqueCounterBits;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

@Deprecated
public class MatchGoalOld extends MultiFileGoal {
	private final ObjectGoal<Database, GSProject> storeGoal;
	protected final ExecutionContext bundle;
	private final Map<String, MatchingResult> matchResults;
	private final Map<File, String> fileToKeyMap;

	private FastqKMerMatcher matcher;
	private Database database;
	private ResultReporter reporter;
	private KMerUniqueCounter uniqueCounter;

	@SafeVarargs
	public MatchGoalOld(GSProject project, GoalKey key, ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapGoal,
						ObjectGoal<Database, GSProject> storeGoal, ExecutionContext bundle, Goal<GSProject>... deps) {
		super(project, key, fastqMapGoal, Goal.append(deps, storeGoal));
		this.storeGoal = storeGoal;
		this.bundle = bundle;
		fileToKeyMap = new HashMap<File, String>();
		matchResults = new HashMap<String, MatchingResult>();
	}

	@Override
	protected void enterFileAndKey(File file, String key) {
		fileToKeyMap.put(file, key);
	}

	@Override
	protected FileType getFileType() {
		return FileType.CSV;
	}

	@Override
	protected void makeFile(File file) {
		try {
			File filteredFile = null;
			File krakenOutStyleFile = null;
			StreamingResourceStream fastqs = fileToFastqs.get(file);

			if (booleanConfigValue(GSConfigKey.WRITED_FILTERED_FASTQ)) {
				filteredFile = getProject().getOutputFile(null, null, file.getName(), FileType.FASTQ_RES,
						true);
			}
			if (booleanConfigValue(GSConfigKey.WRITED_KRAKEN_STYLE_OUT)) {
				krakenOutStyleFile = getProject().getOutputFile(null, null, file.getName(),
						FileType.KRAKEN_OUT_RES, false);
			}

			if (matcher == null) {
				database = storeGoal.get();
				database.getKmerStore().setUseFilter(booleanConfigValue(GSConfigKey.USE_BLOOM_FILTER_FOR_MATCH));

				SmallTaxTree taxTree = database.getTaxTree();
				matcher = createMatcher(database.convertKMerStore(),
						(booleanConfigValue(GSConfigKey.CLASSIFY_READS) && !GSGoalKey.MATCHLR.equals(getKey()))
								? taxTree
								: null,
						bundle, booleanConfigValue(GSConfigKey.WITH_PROBS));
				reporter = new ResultReporter();
				uniqueCounter = booleanConfigValue(GSConfigKey.COUNT_UNIQUE_KMERS)
						? new KMerUniqueCounterBits(database.getKmerStore(),
								intConfigValue(GSConfigKey.MAX_KMER_RES_COUNTS) > 0)
						: null;
			}
			if (uniqueCounter != null) {
				uniqueCounter.clear();
			}
			MatchingResult res = matcher.runMatcher(fastqs, filteredFile, krakenOutStyleFile, uniqueCounter);
			res.completeResults(database);
			storeResult(file, res);
			writeOutputFile(file, res);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected void storeResult(File file, MatchingResult res) {
		matchResults.put(fileToKeyMap.get(file), res);
	}

	public Map<String, MatchingResult> getMatchResults() {
		return matchResults;
	}

	@Override
	public void dump() {
		super.dump();
		dumpMatcher();
	}

	protected void dumpMatcher() {
		if (matcher != null) {
			matcher.dump();
			matcher = null;
		}
	}

	protected void writeOutputFile(File file, MatchingResult result) throws IOException {
		try (PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file))) {
			reporter.printMatchResult(result, out);
		}
	}

	protected FastqKMerMatcher createMatcher(KMerSortedArray<SmallTaxIdNode> store, SmallTaxTree taxTree,
			ExecutionContext bundle, boolean withProbs) {
		return new FastqKMerMatcher(store, intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES),
				intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE), bundle, withProbs, intConfigValue(GSConfigKey.MAX_KMER_RES_COUNTS),
				taxTree, intConfigValue(GSConfigKey.MAX_CLASSIFICATION_PATHS),
				doubleConfigValue(GSConfigKey.MAX_READ_TAX_ERROR_COUNT),
				doubleConfigValue(GSConfigKey.MAX_READ_CLASS_ERROR_COUNT),
				booleanConfigValue(GSConfigKey.KRAKEN_STYLE_MATCH));
	}

	@Override
	protected void endMake() {
		if (matcher != null) {
			matcher.dump();
			matcher = null;
			database = null;
			uniqueCounter = null;
			reporter = null;
		}
	}
}
