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

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.GSFileType;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.FastqKMerMatcher;
import org.metagene.genestrip.match.InlinedFastqKMerMatcher;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerUniqueCounter;
import org.metagene.genestrip.store.KMerUniqueCounterBits;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class MatchResultGoal<P extends GSProject> extends ObjectGoal<Map<String, MatchingResult>, P> {
	private final ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal;
	private final ObjectGoal<Database, P> storeGoal;
	private final ExecutionContext bundle;
	private AfterMatchCallback afterMatchCallback;

	@SafeVarargs
	public MatchResultGoal(P project, GoalKey key, ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal,
                           ObjectGoal<Database, P> storeGoal, ExecutionContext bundle, Goal<P>... deps) {
		super(project, key, Goal.append(deps, fastqMapGoal, storeGoal));
		this.fastqMapGoal = fastqMapGoal;
		this.storeGoal = storeGoal;
		this.bundle = bundle;
	}

	public ObjectGoal<Map<String, StreamingResourceStream>, P> getFastqMapGoal() {
		return fastqMapGoal;
	}

	@Override
	protected void doMakeThis() {
		FastqKMerMatcher matcher = null;
		try {
			Map<String, MatchingResult> matchResults = new HashMap<>();
			Map<String, StreamingResourceStream> map = fastqMapGoal.get();
			Database database = null;
			KMerUniqueCounter uniqueCounter = null;

			for (String key : map.keySet()) {
				File filteredFile = null;
				File krakenOutStyleFile = null;
				StreamingResourceStream fastqs = map.get(key);

				if (booleanConfigValue(GSConfigKey.WRITE_FILTERED_FASTQ)) {
					filteredFile = getProject().getOutputFile(getKey().getName(), key, null, GSFileType.FASTQ_RES,
							true);
				}
				if (booleanConfigValue(GSConfigKey.WRITE_KRAKEN_STYLE_OUT)) {
					krakenOutStyleFile = getProject().getOutputFile(getKey().getName(), key, null,
							GSFileType.KRAKEN_OUT_RES, false);
				}

				if (matcher == null) {
					database = storeGoal.get();
					database.getKmerStore().setUseFilter(booleanConfigValue(GSConfigKey.USE_BLOOM_FILTER_FOR_MATCH));

					SmallTaxTree taxTree = database.getTaxTree();
					matcher = createMatcher(database.convertKMerStore(),
							(booleanConfigValue(GSConfigKey.CLASSIFY_READS) && !GSGoalKey.MATCHRESLR.equals(getKey()))
									? taxTree
									: null,
							bundle, booleanConfigValue(GSConfigKey.WITH_PROBS));
					if (afterMatchCallback != null) {
						matcher.setAfterMatchCallback(new FastqKMerMatcher.AfterMatchCallback() {
							@Override
							public void afterMatch(FastqKMerMatcher.MatcherReadEntry entry, boolean found) {
								afterMatchCallback.afterMatch(entry, found);
							}
						});
					}
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
				matchResults.put(key, res);
				if (afterMatchCallback != null) {
					afterMatchCallback.afterKey(key, res);
				}
			}
			set(matchResults);
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			if (matcher != null) {
				matcher.dump();
			}
		}
	}

	protected FastqKMerMatcher createMatcher(KMerSortedArray<SmallTaxIdNode> store, SmallTaxTree taxTree,
			ExecutionContext bundle, boolean withProbs) {
		if (booleanConfigValue(GSConfigKey.USE_INLINED)) {
			return new InlinedFastqKMerMatcher(store, intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES),
					intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE), bundle, withProbs, intConfigValue(GSConfigKey.MAX_KMER_RES_COUNTS),
					taxTree, intConfigValue(GSConfigKey.MAX_CLASSIFICATION_PATHS),
					doubleConfigValue(GSConfigKey.MAX_READ_TAX_ERROR_COUNT),
					doubleConfigValue(GSConfigKey.MAX_READ_CLASS_ERROR_COUNT),
					booleanConfigValue(GSConfigKey.WRITE_ALL)) {
				@Override
				protected boolean isProgressBar() {
					return booleanConfigValue(GSConfigKey.PROGRESS_BAR);
				}

				@Override
				protected String getProgressBarTaskName() {
					return getKey().getName();
				}			};
		}
		else {
			return new FastqKMerMatcher(store, intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES),
					intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE), bundle, withProbs, intConfigValue(GSConfigKey.MAX_KMER_RES_COUNTS),
					taxTree, intConfigValue(GSConfigKey.MAX_CLASSIFICATION_PATHS),
					doubleConfigValue(GSConfigKey.MAX_READ_TAX_ERROR_COUNT),
					doubleConfigValue(GSConfigKey.MAX_READ_CLASS_ERROR_COUNT),
					booleanConfigValue(GSConfigKey.WRITE_ALL)) {
				@Override
				protected boolean isProgressBar() {
					return booleanConfigValue(GSConfigKey.PROGRESS_BAR);
				}

				@Override
				protected String getProgressBarTaskName() {
					return getKey().getName();
				}
			};
		}
	}

	public void setAfterMatchCallback(AfterMatchCallback afterMatchCallback) {
		this.afterMatchCallback = afterMatchCallback;
	}

	public interface AfterMatchCallback {
		public void afterKey(String key, MatchingResult res);
		public void afterMatch(FastqKMerMatcher.MatcherReadEntry entry, boolean found);
	}
}
