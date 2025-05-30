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
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerSortedArray.UpdateValueProvider;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class DBGoal extends FastaReaderGoal<Database> {
	private final ObjectGoal<AccessionMap, GSProject> accessionTrieGoal;
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final ObjectGoal<Database, GSProject> filledStoreGoal;
	private final boolean minUpdate;

	private KMerSortedArray<String> store;

	@SafeVarargs
	public DBGoal(GSProject project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, GSProject> categoriesGoal,
				  ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
				  ObjectGoal<TaxTree, GSProject> taxTreeGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
				  ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalGoal,
			ObjectGoal<AccessionMap, GSProject> accessionTrieGoal, ObjectGoal<Database, GSProject> filledStoreGoal,
			Goal<GSProject>... deps) {
		super(project, GSGoalKey.UPDATE_DB, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, Goal.append(deps, taxTreeGoal, accessionTrieGoal, filledStoreGoal));
		this.taxTreeGoal = taxTreeGoal;
		this.accessionTrieGoal = accessionTrieGoal;
		this.filledStoreGoal = filledStoreGoal;
		minUpdate = project.booleanConfigValue(GSConfigKey.MIN_UPDATE);
	}
	
	@Override
	protected void doMakeThis() {       
		try {
			Database wrapper = filledStoreGoal.get();
			store = wrapper.getKmerStore();
			readFastas();
			store.fix();
			set(new Database(store, wrapper.getTaxTree()));
			if (getLogger().isTraceEnabled()) {
				getLogger().trace("KMers moved: " + store.getKMersMoved());
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			store = null;
			cleanUpThreads();
		}
	}

	protected AbstractRefSeqFastaReader createFastaReader() {
		return new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES), taxTreeGoal.get(), taxNodesGoal.get(),
				accessionTrieGoal.get(), store, intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
				(Rank) configValue(GSConfigKey.MAX_GENOMES_PER_TAXID_RANK),
				longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID),
				intConfigValue(GSConfigKey.MAX_DUST),
				intConfigValue(GSConfigKey.STEP_SIZE),
				booleanConfigValue(GSConfigKey.UPDATE_WITH_COMPLETE_GENOMES_ONLY));
	}

	protected class MyFastaReader extends AbstractStoreFastaReader {
		private final KMerSortedArray<String> store;
		private final UpdateValueProvider<String> provider;

		public MyFastaReader(int bufferSize, TaxTree taxTree, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, KMerSortedArray<String> store,
							 int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly) {
			super(bufferSize, taxNodes, accessionMap, store.getK(), maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, stepSize, completeGenomesOnly);
			this.store = store;
			provider = new UpdateValueProvider<String>() {
				// Caches for last results of getLeastCommonAncestor()
				private String lastOldValue;
				private TaxIdNode lastNode;
				private String lastLCA;

				@Override
				public String getUpdateValue(String oldValue) {
					// Minimal result cache to improve speed - works 95% of the time.
					if (oldValue == lastOldValue && node == lastNode) {
						return lastLCA;
					}

					TaxIdNode oldNode = taxTree.getNodeByTaxId(oldValue);
					TaxIdNode lcaNode = taxTree.getLeastCommonAncestor(oldNode, node);

					lastOldValue = oldValue;
					lastNode = node;
					lastLCA = lcaNode != null ? lcaNode.getTaxId() : oldValue;

					return lastLCA;
				}
			};
		}

		@Override
		protected void infoLine() {
			if (!ignoreMap) {
				updateNodeFromInfoLine();
			}
			if (minUpdate) {
				// This means we use all regions that overlap deal with our taxids.
				// It might be more than what is in the DB but still less than the entire RefSeq.
				if (node != null && (taxNodes.isEmpty() || taxNodes.contains(node))) {
					includeRegion = true;
				}
			}
			else {
				includeRegion = true;
			}
		}

		@Override
		public boolean isAllowMoreKmers() {
			return true;
		}

		@Override
		protected boolean handleStore() {
			return store.update(byteRingBuffer.getStandardKMer(), provider);
		}
	}
}