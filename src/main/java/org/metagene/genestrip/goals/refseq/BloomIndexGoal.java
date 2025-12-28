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
import org.metagene.genestrip.bloom.AbstractKMerBloomFilter;
import org.metagene.genestrip.bloom.MurmurKMerBloomFilter;
import org.metagene.genestrip.bloom.XORKMerBloomFilter;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerSortedArray.KMerSortedArrayVisitor;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

public class BloomIndexGoal<P extends GSProject> extends ObjectGoal<AbstractKMerBloomFilter, P> {
	private final ObjectGoal<Database, P> filledStoreGoal;

	@SafeVarargs
	public BloomIndexGoal(P project, ObjectGoal<Database, P> filledStoreGoal, Goal<P>... deps) {
		super(project, GSGoalKey.FILL_INDEX, append(deps, filledStoreGoal));
		this.filledStoreGoal = filledStoreGoal;
	}

	@Override
	protected void doMakeThis() {
		Database wrapper = filledStoreGoal.get();
		KMerSortedArray<String> store = wrapper.getKmerStore();
		SmallTaxTree taxTree = wrapper.getTaxTree();

		long[] counter = new long[1];
		store.visit(new KMerSortedArrayVisitor<String>() {
			@Override
			public void nextValue(KMerSortedArray<String> trie, long kmer, int index, long i) {
				String taxid = store.getValueForIndex(index);
				if (taxid != null) {
					SmallTaxIdNode node = taxTree.getNodeByTaxId(taxid);
					if (node.isRequested()) {
						counter[0]++;
					}
				}
			}
		});

		if (getLogger().isInfoEnabled()) {
			getLogger().info("Bloom filter size: " + counter[0] + " entries.");
		}

		AbstractKMerBloomFilter filter = booleanConfigValue(GSConfigKey.XOR_BLOOM_HASH) ?
				new XORKMerBloomFilter(intConfigValue(GSConfigKey.KMER_SIZE),
						doubleConfigValue(GSConfigKey.BLOOM_FILTER_FPP)) :
				new MurmurKMerBloomFilter(intConfigValue(GSConfigKey.KMER_SIZE),
				doubleConfigValue(GSConfigKey.BLOOM_FILTER_FPP));
		filter.ensureExpectedSize(counter[0], false);

		store.visit(new KMerSortedArrayVisitor<String>() {
			@Override
			public void nextValue(KMerSortedArray<String> trie, long kmer, int index, long i) {
				String taxid = store.getValueForIndex(index);
				if (taxid != null) {
					SmallTaxIdNode node = taxTree.getNodeByTaxId(taxid);
					if (node.isRequested()) {
						filter.putLong(kmer);
					}
				}
			}
		});
		set(filter);
	}
}