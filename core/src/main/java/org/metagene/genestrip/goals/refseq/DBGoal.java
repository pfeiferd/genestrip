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
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
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
import org.metagene.genestrip.refseq.ReworkingStoreFastaReader;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.store.KMerStore.UpdateValueProvider;
import org.metagene.genestrip.store.RadixKMerStore;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * Goal ({@code UPDATE_DB}) that updates an already-filled database by re-reading the RefSeq (and
 * additional) FASTA files and merging each k-mer's stored taxid with the k-mer's new node via their
 * lowest common ancestor, then stamps Genestrip version, title and creation-date properties.
 * Produces the updated {@link Database}.
 *
 * @param <P> the project type
 */
public class DBGoal<P extends GSProject> extends FastaReaderGoal<Database, P> {
	private final ObjectGoal<AccessionMap, P> accessionMapGoal;
	private final ObjectGoal<TaxTree, P> taxTreeGoal;
	private final ObjectGoal<Database, P> filledStoreGoal;
	private final boolean minUpdate;

	private KMerStore<String> store;

	/**
	 * Creates the goal, wiring the tax tree, accession map and filled-store goals it updates in place.
	 *
	 * @param project the project this goal belongs to
	 * @param bundle the execution context providing worker threads
	 * @param categoriesGoal the goal supplying the RefSeq categories to include
	 * @param taxNodesGoal the goal supplying the set of taxonomy nodes
	 * @param taxTreeGoal the goal supplying the taxonomy tree
	 * @param fnaFilesGoal the goal supplying the downloaded RefSeq FASTA files
	 * @param additionalGoal the goal supplying additional FASTA files mapped to taxonomy nodes
	 * @param accessionMapGoal the goal supplying the accession-to-taxid map
	 * @param filledStoreGoal the goal supplying the filled database to update
	 * @param deps additional goals this goal depends on
	 */
	@SafeVarargs
	public DBGoal(P project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
				  ObjectGoal<Set<TaxIdNode>, P> taxNodesGoal,
				  ObjectGoal<TaxTree, P> taxTreeGoal, RefSeqFnaFilesDownloadGoal<P> fnaFilesGoal,
				  ObjectGoal<Map<File, TaxIdNode>, P> additionalGoal,
			ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<Database, P> filledStoreGoal,
			Goal<P>... deps) {
		super(project, GSGoalKey.UPDATE_DB, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, Goal.append(deps, taxTreeGoal, accessionMapGoal, filledStoreGoal));
		this.taxTreeGoal = taxTreeGoal;
		this.accessionMapGoal = accessionMapGoal;
		this.filledStoreGoal = filledStoreGoal;
		minUpdate = project.booleanConfigValue(GSConfigKey.MIN_UPDATE);
	}

	@Override
	protected boolean isIncludeRefSeqFna() {
		return true;
	}

	@Override
	protected void doMakeThis() {       
		try {
			Database wrapper = filledStoreGoal.get();
			store = wrapper.getKmerStore();
			readFastas();
			// readFastas() reassigned k-mer values via the bulk update paths, which do not touch the
			// per-taxid count cache; drop it so it is recomputed on the next read (and baked into the
			// serialized database, see AbstractKMerStore#writeObject).
			store.invalidateNKmersPerTaxid();
			String gsVersion = GSProject.getGenestripRuntimeVersion();
			if (gsVersion != null) {
				getProject().setAdditionalProperty(GSProject.GENESTRIP_VERSION, gsVersion);
			}
			String gsTitle = GSProject.getGenestripRuntimeTitle();
			if (gsTitle != null) {
				getProject().setAdditionalProperty(GSProject.GENESTRIP_TITLE, gsTitle);
			}
			DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
			getProject().setAdditionalProperty(GSProject.DB_CREATION_DATE, dateFormat.format(new Date()));
			set(new Database(store, wrapper.getTaxTree(), getProject().getAllAsProperties()));
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			store = null;
			cleanUpThreads();
		}
	}

	protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
		// Lookup mode (createNodes = false): the artificial data/file/id nodes were created during the
		// fill; the update only looks them up (using the same key computation as their creation).
		return new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES), taxTreeGoal.get(), taxNodesGoal.get(),
				accessionMapGoal.get(), store, intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
				(Rank) configValue(GSConfigKey.MAX_GENOMES_PER_TAXID_RANK),
				longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID),
				intConfigValue(GSConfigKey.MAX_DUST),
				intConfigValue(GSConfigKey.STEP_SIZE),
				booleanConfigValue(GSConfigKey.UPDATE_WITH_COMPLETE_GENOMES_ONLY),
				null,
				booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES),
				booleanConfigValue(GSConfigKey.DATA_NODES),
				booleanConfigValue(GSConfigKey.FILE_NODES),
				booleanConfigValue(GSConfigKey.ID_NODES));
	}

	/**
	 * FASTA reader that updates each existing k-mer's taxid to the lowest common ancestor of its
	 * current node and the k-mer's new node.
	 */
	protected class MyFastaReader extends ReworkingStoreFastaReader {
		// Number of k-mers gathered before a batched flush; sized to overlap enough independent
		// cache-missing lookups without the per-batch bookkeeping outweighing the memory-level
		// parallelism it buys (see RadixKMerStore.updateBatch).
		private static final int UPDATE_BATCH_SIZE = 128;

		private final KMerStore<String> store;
		// Non-null exactly when the store is a RadixKMerStore, in which case k-mers are updated in
		// batches to expose memory-level parallelism instead of one at a time.
		private final RadixKMerStore<String> radixStore;
		private final RadixKMerStore.BatchBuffers batch;
		private final UpdateValueProvider<String> provider;

		/**
		 * Creates the reader that updates existing k-mers to the lowest common ancestor of their
		 * current and new nodes.
		 *
		 * @param bufferSize the FASTA read buffer size in bytes
		 * @param taxTree the taxonomy tree
		 * @param taxNodes the set of taxonomy nodes to consider
		 * @param accessionMap the accession-to-taxid map
		 * @param store the k-mer store to update
		 * @param maxGenomesPerTaxId the maximum number of genomes per taxid
		 * @param maxGenomesPerTaxIdRank the rank at which the genome limit applies
		 * @param maxKmersPerTaxId the maximum number of k-mers per taxid
		 * @param maxDust the maximum dust (low-complexity) threshold
		 * @param stepSize the k-mer sampling step size
		 * @param completeGenomesOnly whether to restrict to complete genomes only
		 * @param regionsPerTaxid the trie of regions per taxid
		 * @param enableLowerCaseBases whether lowercase bases are treated as valid
		 * @param dataNodes whether artificial {@code DATA} nodes are used
		 * @param fileNodes whether artificial {@code FILE} nodes are used
		 * @param idNodes whether artificial {@code ID} nodes are used
		 */
		@SuppressWarnings("unchecked")
		public MyFastaReader(int bufferSize, TaxTree taxTree, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, KMerStore<String> store,
							 int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases,
							 boolean dataNodes, boolean fileNodes, boolean idNodes) {
			// Lookup mode: the fill already created every artificial node, so no id generator is needed.
			super(bufferSize, taxNodes, accessionMap, store.getK(), maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, stepSize, completeGenomesOnly, regionsPerTaxid, enableLowerCaseBases,
					taxTree, dataNodes, fileNodes, idNodes, false, null);
			this.store = store;
			if (store instanceof RadixKMerStore) {
				radixStore = (RadixKMerStore<String>) store;
				batch = new RadixKMerStore.BatchBuffers(UPDATE_BATCH_SIZE);
			} else {
				radixStore = null;
				batch = null;
			}
			provider = new UpdateValueProvider<>() {
				// Caches for last results of getLowestCommonAncestor()
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
					TaxIdNode lcaNode = taxTree.getLowestCommonAncestor(oldNode, node);

					lastOldValue = oldValue;
					lastNode = node;
					lastLCA = lcaNode != null ? lcaNode.getTaxId() : oldValue;

					return lastLCA;
				}
			};
		}

		@Override
		protected void infoLine() {
			if (ignoreMap) {
				node = mappedNode;
			}
			else {
				updateNodeFromInfoLine();
			}

			if (minUpdate) {
				// This means we use all regions that overlap deal with our taxids.
				// It might be more than what is in the DB but still less than the entire RefSeq.
				if (node != null && (taxNodes.isEmpty() || taxNodes.contains(node))) {
					includeRegion = true;
					node = reworkNode();
				}
			}
			else {
				includeRegion = true;
				if (node != null) {
					node = reworkNode();
				}
			}
		}

		@Override
		protected void endRegion() {
			// Flush the current region's pending k-mers before infoLine() moves 'node' to the next
			// region: every batch must consist of k-mers that share the region's node, since the
			// provider merges against that node. A region is the coarsest safe flush boundary; the
			// buffer is also flushed mid-region once it fills up (see handleStore()).
			if (batch != null && !batch.isEmpty()) {
				radixStore.updateBatch(batch, provider);
			}
		}

		@Override
		public boolean isAllowMoreKmers() {
			return true;
		}

		@Override
		protected boolean handleStore() {
			if (batch != null) {
				if (batch.add(byteRingBuffer.getStandardKMer())) {
					radixStore.updateBatch(batch, provider);
				}
				// The counted-k-mer return value is unused by this reader (endRegion does no region
				// bookkeeping and isAllowMoreKmers() is always true), so the deferred move result of a
				// batched k-mer need not be reported here.
				return false;
			}
			return store.update(byteRingBuffer.getStandardKMer(), provider);
		}
	}
}