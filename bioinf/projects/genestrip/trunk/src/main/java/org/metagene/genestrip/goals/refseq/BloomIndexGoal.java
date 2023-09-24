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
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerSortedArray.KMerSortedArrayVisitor;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.store.KMerStoreWrapper;

public class BloomIndexGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final ObjectGoal<KMerStoreWrapper, GSProject> filledStoreGoal;
	private BloomIndexedGoal indexedGoal;

	@SafeVarargs
	public BloomIndexGoal(GSProject project, String name, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, ObjectGoal<KMerStoreWrapper, GSProject> filledStoreGoal,
			Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.SER),
				append(deps, taxTreeGoal, taxNodesGoal, filledStoreGoal));
		this.taxTreeGoal = taxTreeGoal;
		this.taxNodesGoal = taxNodesGoal;
		this.filledStoreGoal = filledStoreGoal;
	}
	
	public void setBloomIndexedGoal(BloomIndexedGoal indexedGoal) {
		this.indexedGoal = indexedGoal;
	}
	

	@Override
	public void makeFile(File filterFile) {
		try {
			KMerStoreWrapper wrapper = filledStoreGoal.get();
			KMerSortedArray<String> store = wrapper.getKmerStore();
			TaxTree taxTree = taxTreeGoal.get();
			Set<TaxIdNode> nodes = taxNodesGoal.get();

			long[] counter = new long[1];
			store.visit(new KMerSortedArrayVisitor<String>() {
				@Override
				public void nextValue(KMerSortedArray<String> trie, long kmer, short index, long i) {
					String taxid = store.getValueForIndex(index);
					if (taxid != null) {
						TaxIdNode node = taxTree.getNodeByTaxId(taxid);
						if (nodes.contains(node)) {
							counter[0]++;
						}
					}
				}
			});

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Bloom filter size " + counter[0]);
			}

			MurmurCGATBloomFilter filter = new MurmurCGATBloomFilter(getProject().getConfig().getKMerSize(),
					getProject().getConfig().getKMerFastBloomFpp());
			filter.ensureExpectedSize(counter[0], false);

			
			store.visit(new KMerSortedArrayVisitor<String>() {
				@Override
				public void nextValue(KMerSortedArray<String> trie, long kmer, short index, long i) {
					String taxid = store.getValueForIndex(index);
					if (taxid != null) {
						TaxIdNode node = taxTree.getNodeByTaxId(taxid);
						if (nodes.contains(node)) {
							filter.putLong(kmer);
						}
					}
				}
			});
			
			filter.save(filterFile);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File saved " + filterFile);
			}
			
			indexedGoal.setFilter(filter);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}