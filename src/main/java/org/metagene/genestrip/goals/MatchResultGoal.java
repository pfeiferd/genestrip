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
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.FastqKMerMatcher;
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

public class MatchResultGoal extends ObjectGoal<Map<String, MatchingResult>, GSProject> {
	private final ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapGoal;
	private final ObjectGoal<Database, GSProject> storeGoal;
	protected final ExecutionContext bundle;

	@SafeVarargs
	public MatchResultGoal(GSProject project, GoalKey key, ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapGoal,
                           ObjectGoal<Database, GSProject> storeGoal, ExecutionContext bundle, Goal<GSProject>... deps) {
		super(project, key, Goal.append(deps, fastqMapGoal, storeGoal));
		this.fastqMapGoal = fastqMapGoal;
		this.storeGoal = storeGoal;
		this.bundle = bundle;
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

				if (booleanConfigValue(GSConfigKey.WRITED_FILTERED_FASTQ)) {
					filteredFile = getProject().getOutputFile(getKey().getName(), key, null, FileType.FASTQ_RES,
							true);
				}
				if (booleanConfigValue(GSConfigKey.WRITED_KRAKEN_STYLE_OUT)) {
					krakenOutStyleFile = getProject().getOutputFile(getKey().getName(), key, null,
							FileType.KRAKEN_OUT_RES, false);
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
		return new FastqKMerMatcher(store, intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES),
				intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE), bundle, withProbs, intConfigValue(GSConfigKey.MAX_KMER_RES_COUNTS),
				taxTree, intConfigValue(GSConfigKey.MAX_CLASSIFICATION_PATHS),
				doubleConfigValue(GSConfigKey.MAX_READ_TAX_ERROR_COUNT),
				doubleConfigValue(GSConfigKey.MAX_READ_CLASS_ERROR_COUNT));
	}
}
