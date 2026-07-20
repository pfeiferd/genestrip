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

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.bloom.BlockedKMerBloomFilter;
import org.metagene.genestrip.bloom.KMerProbFilter;
import org.metagene.genestrip.bloom.MurmurKMerBloomFilter;
import org.metagene.genestrip.bloom.XORKMerBloomFilter;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.store.KMerStore.IndexedKMerStoreVisitor;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

/**
 * Goal that builds the index Bloom filter ({@link KMerProbFilter}) from the filled database: it
 * inserts every k-mer whose stored taxid belongs to a requested tax node, sizing the filter from
 * the configured false-positive probability.
 *
 * @param <P> the project type
 */
public class BloomIndexGoal<P extends GSProject> extends ObjectGoal<KMerProbFilter, P> {
	private final ObjectGoal<Database, P> filledStoreGoal;

	/**
	 * Creates the goal, wiring the filled database goal whose k-mers are indexed.
	 *
	 * @param project         the project this goal belongs to
	 * @param filledStoreGoal the goal supplying the filled database whose k-mers are indexed
	 * @param deps            any further goals this goal depends on
	 */
	@SafeVarargs
	public BloomIndexGoal(P project, ObjectGoal<Database, P> filledStoreGoal, Goal<P>... deps) {
		super(project, GSGoalKey.FILL_INDEX, append(deps, filledStoreGoal));
		this.filledStoreGoal = filledStoreGoal;
	}

	@Override
	protected void doMakeThis() {
		Database wrapper = filledStoreGoal.get();
		KMerStore<String> store = wrapper.getKmerStore();
		SmallTaxTree taxTree = wrapper.getTaxTree();

		long[] counter = new long[1];
		store.visit(new IndexedKMerStoreVisitor<String>() {
			@Override
			public void nextValue(KMerStore<String> trie, long kmer, int index, long i) {
				String taxid = store.getValueForIndex(index);
				if (taxid != null) {
					SmallTaxIdNode node = taxTree.getNodeByTaxId(taxid);
					if (node != null && node.isRequested()) {
						counter[0]++;
					}
				}
			}
		});

		if (getLogger().isInfoEnabled()) {
			getLogger().info("Bloom filter size: " + counter[0] + " entries.");
		}

		double fpp = doubleConfigValue(GSConfigKey.INDEX_BLOOM_FILTER_FPP);
		KMerProbFilter filter;
		if (fpp == BlockedKMerBloomFilter.DEFAULT_FPP) {
			filter = new BlockedKMerBloomFilter(counter[0]);
		}
		else {
			filter = booleanConfigValue(GSConfigKey.XOR_BLOOM_HASH) ?
					new XORKMerBloomFilter(fpp, counter[0]) : new MurmurKMerBloomFilter(fpp, counter[0]);
		}

		store.visit(new IndexedKMerStoreVisitor<String>() {
			@Override
			public void nextValue(KMerStore<String> trie, long kmer, int index, long i) {
				String taxid = store.getValueForIndex(index);
				if (taxid != null) {
					SmallTaxIdNode node = taxTree.getNodeByTaxId(taxid);
					if (node != null && node.isRequested()) {
						filter.putLong(kmer);
					}
				}
			}
		});
		set(filter);
	}
}