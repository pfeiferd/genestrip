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
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class FillDBGoal extends FastaReaderGoal<Database> {
	private final ObjectGoal<AccessionMap, GSProject> accessionMapGoal;
	private final ObjectGoal<MurmurCGATBloomFilter, GSProject> bloomFilterGoal;
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;

	@SafeVarargs
	public FillDBGoal(GSProject project, ObjectGoal<Set<RefSeqCategory>, GSProject> categoriesGoal,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
					  ObjectGoal<TaxTree, GSProject> taxTreeGoal,
					  RefSeqFnaFilesDownloadGoal fnaFilesGoal,
			ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalGoal,
			ObjectGoal<AccessionMap, GSProject> accessionMapGoal,
			ObjectGoal<MurmurCGATBloomFilter, GSProject> bloomFilterGoal,
			Goal<GSProject>... deps) {
		super(project, GSGoalKey.FILL_DB, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, Goal.append(deps, bloomFilterGoal));
		this.accessionMapGoal = accessionMapGoal;
		this.bloomFilterGoal = bloomFilterGoal;
		this.taxTreeGoal = taxTreeGoal;
	}

	@Override
	protected void doMakeThis() {
		KMerSortedArray<String> store = new KMerSortedArray<>(intConfigValue(GSConfigKey.KMER_SIZE),
				doubleConfigValue(GSConfigKey.BLOOM_FILTER_FPP), null, false);
		// We have to account for the missing entries in the bloom filter due to
		// inherent FPP.
		// This works really well, so we can allow for a low FPP for the bloom filter
		// itself and save memory during db construction.
		// It is a very conservative estimate too, since collisions occur in the process
		// of filling (as opposed to the FPP formula considers a filled
		// bloom filter).
		store.initSize((long) (bloomFilterGoal.get().getEntries()
				* (1 + doubleConfigValue(GSConfigKey.TEMP_BLOOM_FILTER_FPP))));

		try {
			boolean refSeqDB = booleanConfigValue(GSConfigKey.REF_SEQ_DB);

			MyFastaReader fastaReader = new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
					taxNodesGoal.get(), refSeqDB ? accessionMapGoal.get() : null, store,
					intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
					(Rank) configValue(GSConfigKey.MAX_GENOMES_PER_TAXID_RANK),
					longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID),
					intConfigValue(GSConfigKey.MAX_DUST),
					intConfigValue(GSConfigKey.STEP_SIZE));
			readFastas(fastaReader);
			if (getLogger().isWarnEnabled() && fastaReader.tooManyCounter > 0) {
				getLogger().warn("Not stored kmers: " + fastaReader.tooManyCounter);
			}
			TaxTree taxTree = taxTreeGoal.get();
			Iterator<String> taxIt = store.getValues();
			while (taxIt.hasNext()) {
				TaxIdNode node = taxTree.getNodeByTaxId(taxIt.next());
				if (node != null) {
					node.markRequired();
				}
			}
			SmallTaxTree smallTaxTree = taxTreeGoal.get().toSmallTaxTree();
			for (TaxIdNode node : taxNodesGoal.get()) {
				SmallTaxIdNode smallNode = smallTaxTree.getNodeByTaxId(node.getTaxId());
				if (smallNode != null) {
					smallNode.setRequested(true);
					// smallNode.setStoreIndex(store.getIndexForValue(node.getTaxId()));
				}
			}
			ensureAllTreeNodesInDB(smallTaxTree.getRoot(), store);
			store.optimize();
			Database wrapper = new Database(store, smallTaxTree);
			set(wrapper);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected void ensureAllTreeNodesInDB(SmallTaxIdNode node, KMerSortedArray<String> store) {
		if (node == null || node.getTaxId() == null) {
			return;
		}
		store.getAddValueIndex(node.getTaxId());
		SmallTaxIdNode[] subnodes = node.getSubNodes();
		if (subnodes != null) {
			for (int i = 0; i < subnodes.length; i++) {
				ensureAllTreeNodesInDB(subnodes[i], store);
			}
		}
	}

	protected static class MyFastaReader extends AbstractStoreFastaReader {
		private final KMerSortedArray<String> store;
		private long tooManyCounter;

		public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap,
							 KMerSortedArray<String> store, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize) {
			super(bufferSize, taxNodes, accessionMap, store.getK(), maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, stepSize);
			this.store = store;
		}

		@Override
		protected boolean handleStore() {
			if (store.isFull()) {
				tooManyCounter++;
			} else {
				if (node.getTaxId() != null) {
					return store.putLong(byteRingBuffer.getKMer(), byteRingBuffer.getReverseKMer(), node.getTaxId());
				} else {
					if (getLogger().isWarnEnabled()) {
						getLogger().warn("Tax id node without taxid: " + node.getName());
					}
				}
			}
			return false;
		}
	}
}