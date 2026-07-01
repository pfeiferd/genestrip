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
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.store.KMerUniqueCounterBits;
import org.metagene.genestrip.store.TunableKMerStore;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Runs the k-mer matcher over each input fastq stream against the database and produces the per-key
 * {@link MatchingResult}s, optionally also writing filtered and kraken-style output files.
 *
 * @param <P> the project type
 */
public class MatchResultGoal<P extends GSProject> extends ObjectGoal<Map<String, MatchingResult>, P> {
	private final ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal;
	private final ObjectGoal<Database, P> storeGoal;
	private final ExecutionContext bundle;
	private AfterMatchCallback afterMatchCallback;

	/**
	 * Creates the goal wiring the fastq map and database goals that supply the matcher's inputs.
	 *
	 * @param project the project this goal belongs to
	 * @param key the goal key identifying this goal
	 * @param fastqMapGoal the goal supplying the per-key fastq input streams
	 * @param storeGoal the goal supplying the database to match against
	 * @param bundle the execution context providing threads and progress reporting
	 * @param deps additional goals this goal depends on
	 */
	@SafeVarargs
	public MatchResultGoal(P project, GoalKey key, ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal,
                           ObjectGoal<Database, P> storeGoal, ExecutionContext bundle, Goal<P>... deps) {
		super(project, key, Goal.append(deps, fastqMapGoal, storeGoal));
		this.fastqMapGoal = fastqMapGoal;
		this.storeGoal = storeGoal;
		this.bundle = bundle;
	}

	/**
	 * Returns the goal supplying the per-key fastq input streams.
	 *
	 * @return the fastq map goal
	 */
	public ObjectGoal<Map<String, StreamingResourceStream>, P> getFastqMapGoal() {
		return fastqMapGoal;
	}

	@Override
	protected void doMakeThis() {
		FastqKMerMatcher matcher = null;
		try {
			Map<String, MatchingResult> matchResults = new LinkedHashMap<>(); // To preserver order of keys.
			Map<String, StreamingResourceStream> map = fastqMapGoal.get();
			Database database = null;
			KMerUniqueCounterBits uniqueCounter = null;

			for (String key : map.keySet()) {
				File filteredFile = null;
				File krakenOutStyleFile = null;
				StreamingResourceStream fastqs = map.get(key);

				if (booleanConfigValue(GSConfigKey.WRITE_FILTERED_FASTQ)) {
					filteredFile = getProject().getOutputFile(getKey().getName(), key, null, GSFileType.FASTQ_RES,
							booleanConfigValue(GSConfigKey.GZIP_FASTQ_OUTPUT));
				}
				if (booleanConfigValue(GSConfigKey.WRITE_KRAKEN_STYLE_OUT)) {
					krakenOutStyleFile = getProject().getOutputFile(getKey().getName(), key, null,
							GSFileType.KRAKEN_OUT_RES, false);
				}

				if (matcher == null) {
					database = storeGoal.get();
					SmallTaxTree taxTree = database.getTaxTree();
					KMerStore<SmallTaxIdNode> store = database.convertKMerStore();

					// The pre-filter is optional (see TunableKMerStore); apply it only when supported.
					if (store instanceof TunableKMerStore) {
						((TunableKMerStore<SmallTaxIdNode>) store)
								.setUseFilter(booleanConfigValue(GSConfigKey.USE_BLOOM_FILTER_FOR_MATCH));
					}

					matcher = createMatcher(store,
							(booleanConfigValue(GSConfigKey.CLASSIFY_READS) && !GSGoalKey.MATCHRESLR.equals(getKey()))
									? taxTree
									: null,
							bundle, booleanConfigValue(GSConfigKey.WITH_PROBS),
							database.getConfigInfo().getProperty(GSProject.DB_MD5));
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

	/**
	 * Creates the {@link FastqKMerMatcher} for the given store and tax tree, configured from the goal's config
	 * values.
	 *
	 * @param store the k-mer store to match reads against
	 * @param taxTree the tax tree used for read classification, or {@code null} to disable classification
	 * @param bundle the execution context providing threads and progress reporting
	 * @param withProbs whether match probabilities are computed
	 * @param dbMD5 the MD5 checksum of the database used
	 * @return the configured matcher
	 */
	protected FastqKMerMatcher createMatcher(KMerStore<SmallTaxIdNode> store, SmallTaxTree taxTree,
			ExecutionContext bundle, boolean withProbs, String dbMD5) {
		return new FastqKMerMatcher(store, intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES),
				intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE), bundle, withProbs, intConfigValue(GSConfigKey.MAX_KMER_RES_COUNTS),
				taxTree, intConfigValue(GSConfigKey.MAX_CLASSIFICATION_PATHS),
				doubleConfigValue(GSConfigKey.MAX_READ_TAX_ERROR_COUNT),
				doubleConfigValue(GSConfigKey.MAX_READ_CLASS_ERROR_COUNT),
				booleanConfigValue(GSConfigKey.WRITE_ALL),
				intConfigValue(GSConfigKey.MIN_KMERS_FOR_CLASS),
				dbMD5) {
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

	/**
	 * Sets the callback invoked during matching.
	 *
	 * @param afterMatchCallback the callback to invoke, or {@code null} to disable it
	 */
	public void setAfterMatchCallback(AfterMatchCallback afterMatchCallback) {
		this.afterMatchCallback = afterMatchCallback;
	}

	/**
	 * Callback invoked during matching: once per finished key and once per matched read.
	 */
	public interface AfterMatchCallback {
		/**
		 * Called once all reads for a map key have been matched.
		 *
		 * @param key the map key that was completed
		 * @param res the matching result for that key
		 */
		public void afterKey(String key, MatchingResult res);

		/**
		 * Called for each matched read.
		 *
		 * @param entry the read entry that was matched
		 * @param found whether the read matched the database
		 */
		public void afterMatch(FastqKMerMatcher.MatcherReadEntry entry, boolean found);
	}
}
