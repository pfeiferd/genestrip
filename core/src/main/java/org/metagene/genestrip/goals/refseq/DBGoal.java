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
import org.metagene.genestrip.tax.TaxNodeSelection;

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
	private final GSConfigKey.UpdateScope updateScope;
	/** The tax ids selected for the database, with their descendants; held for the pass. */
	private Set<TaxIdNode> selectedTaxNodes;
	/** The tax ids `taxids.txt' struck out, with their descendants; held for the pass. */
	private Set<TaxIdNode> excludedTaxNodes;

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
				  ObjectGoal<TaxNodeSelection, P> taxNodesGoal,
				  ObjectGoal<TaxTree, P> taxTreeGoal, RefSeqFnaFilesDownloadGoal<P> fnaFilesGoal,
				  ObjectGoal<Map<File, TaxIdNode>, P> additionalGoal,
			ObjectGoal<AccessionMap, P> accessionMapGoal, ObjectGoal<Database, P> filledStoreGoal,
			Goal<P>... deps) {
		// true, and not the configured `refseq.filldb': this pass reads the RefSeq whether or not the
		// database was filled from it. A k-mer's tax id is the lowest common ancestor of the taxa of
		// every genome containing it, and the genomes that would raise it above the requested taxa are
		// precisely the ones of other taxa - which live in the RefSeq release. Skipping it here would
		// leave every k-mer that a relative also carries claimed for the requested taxon, and reads of
		// that relative would be reported as the requested one.
		super(project, GSGoalKey.UPDATE_DB, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, true, Goal.append(deps, taxTreeGoal, accessionMapGoal, filledStoreGoal));
		this.taxTreeGoal = taxTreeGoal;
		this.accessionMapGoal = accessionMapGoal;
		this.filledStoreGoal = filledStoreGoal;
		updateScope = (GSConfigKey.UpdateScope) project.configValue(GSConfigKey.UPDATE_SCOPE);
	}

	/**
	 * Decides whether a contig takes part in the update, given the configured scope, the node the
	 * contig resolved to, and whether it comes from the RefSeq release.
	 * <p>
	 * The scope restricts the <em>release</em> and nothing else, which is what its name
	 * {@code refseq.updateScope} says: a genome the project supplies itself - an entry of
	 * {@code additional.txt} or one downloaded from Genbank - always takes part. Those genomes are
	 * the ones the database was filled from, each under its own identity, so the update meets every
	 * one of them under the identity the fill gave it and {@code LCA(n, n) = n} leaves it where it
	 * is. What it does do there is the half of the update that a per-assembly database still needs:
	 * a k-mer that two of those genomes share settles on their common ancestor instead of staying
	 * claimed for whichever of them happened to be read first. Only the release presents the same
	 * genome a second time under another identity, so only the release is worth restricting.
	 * <p>
	 * A {@code null} node means the contig's accession is not in the accession map, so its taxon is
	 * unknown. Such a contig cannot raise anything - the lowest common ancestor with no node is the
	 * stored value - and it is treated as belonging to no selected taxon, which is what
	 * {@link GSConfigKey.UpdateScope#OTHER_TAXA_ONLY} says about it and why
	 * {@link GSConfigKey.UpdateScope#OWN_TAXA_ONLY} leaves it out.
	 *
	 * @param node the node the contig resolved to, or {@code null} if its taxon is unknown
	 * @param fromRefSeqRelease whether the contig stems from the RefSeq release rather than from a
	 *                          fasta the project supplies itself
	 * @return whether the contig's k-mers are to be merged into the store
	 */
	boolean isContigInScope(TaxIdNode node, boolean fromRefSeqRelease) {
		if (!fromRefSeqRelease) {
			return true;
		}
		switch (updateScope) {
		case OWN_TAXA_ONLY:
			// An empty selection means no restriction, so every known contig qualifies as "ours".
			return node != null && (selectedTaxNodes.isEmpty() || selectedTaxNodes.contains(node));
		case OTHER_TAXA_ONLY:
			// No special case for an empty selection is needed here: nothing is then "ours", so
			// nothing is skipped, which is the same no-restriction reading as above.
			return !selectedTaxNodes.contains(node);
		case ALL_BUT_EXCLUDED:
			// Like ALL, except for the branches `taxids.txt' struck out. A contig whose taxon is
			// unknown is not one of them and takes part, as it does under ALL: it cannot raise
			// anything anyway, the lowest common ancestor with no node being the stored value.
			return !excludedTaxNodes.contains(node);
		default:
			return true;
		}
	}

	@Override
	protected void doMakeThis() {       
		try {
			Database wrapper = filledStoreGoal.get();
			store = wrapper.getKmerStore();
			// Once for the whole pass. isContigInScope() is asked for every contig of the release, and
			// a goal is not a lookup table to be consulted a few hundred thousand times.
			TaxNodeSelection selection = taxNodesGoal.get();
			selectedTaxNodes = selection.getSelected();
			excludedTaxNodes = selection.getExcluded();
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
			selectedTaxNodes = null;
			excludedTaxNodes = null;
		}
	}

	protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid, AbstractRefSeqFastaReader.GenomeKeyTrie admittedGenomes) {
		// Lookup mode (createNodes = false): the artificial data/file/id nodes were created during the
		// fill; the update only looks them up (using the same key computation as their creation).
		// The per-taxon limits maxGenomesPerTaxid, maxPerTaxidRank and maxKMersPerTaxid are
		// deliberately not passed on: this reader overrides infoLine() and isAllowMoreKmers(), the
		// only two places that read them, so they would have no effect. Capping the update would be
		// wrong in any case - a k-mer is raised to the common ancestor of every genome carrying it,
		// and a genome left out of the fill still carries it.
		return new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES), taxTreeGoal.get(), taxNodesGoal.get().getSelected(),
				accessionMapGoal.get(), store,
				intConfigValue(GSConfigKey.MAX_DUST),
				intConfigValue(GSConfigKey.KMER_SAMPLING),
				booleanConfigValue(GSConfigKey.UPDATE_WITH_ASSEMBLY_ACCESSIONS_ONLY),
				null,
				null,
				booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES),
				booleanConfigValue(GSConfigKey.DATA_NODES),
				booleanConfigValue(GSConfigKey.FILE_NODES),
				booleanConfigValue(GSConfigKey.ID_NODES),
				(Rank) configValue(GSConfigKey.FOLD_TAXA_BELOW));
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
		 * @param maxDust the maximum dust (low-complexity) threshold
		 * @param kMerSampling the k-mer sampling step size
		 * @param assemblyAccessionsOnly whether only genomic accessions are considered, dropping `NG_`, `NT_` and `NW_`
		 * @param contigsPerTaxid the trie of contigs per taxid
		 * @param admittedGenomes the shared set of genome keys admitted so far
		 * @param enableLowerCaseBases whether lowercase bases are treated as valid
		 * @param dataNodes whether artificial {@code DATA} nodes are used
		 * @param fileNodes whether artificial {@code FILE} nodes are used
		 * @param idNodes whether artificial {@code ID} nodes are used
		 */
		@SuppressWarnings("unchecked")
		public MyFastaReader(int bufferSize, TaxTree taxTree, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, KMerStore<String> store,
							 int maxDust, int kMerSampling, boolean assemblyAccessionsOnly, StringLong2DigitTrie contigsPerTaxid, AbstractRefSeqFastaReader.GenomeKeyTrie admittedGenomes, boolean enableLowerCaseBases,
							 boolean dataNodes, boolean fileNodes, boolean idNodes, Rank foldTaxaBelow) {
			// Lookup mode: the fill already created every artificial node, so no id generator is needed.
			// The per-taxon genome and k-mer limits are handed to the superclass as their no-limit
			// values: this reader overrides infoLine() and isAllowMoreKmers(), which are the only
			// readers of them, so no limit can take effect here and pretending otherwise would only
			// mislead. See the comment at the call site in createFastaReader().
			super(bufferSize, taxNodes, accessionMap, store.getK(), Integer.MAX_VALUE, null, Long.MAX_VALUE, maxDust, kMerSampling, assemblyAccessionsOnly, contigsPerTaxid, admittedGenomes, enableLowerCaseBases,
					taxTree, dataNodes, fileNodes, idNodes, false, null, foldTaxaBelow);
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

			if (isContigInScope(node, isRefSeqReleaseContig())) {
				includeContig = true;
				if (node != null) {
					node = reworkNode();
				}
			}
		}

		@Override
		protected void endContig() {
			// Flush the current contig's pending k-mers before infoLine() moves 'node' to the next
			// contig: every batch must consist of k-mers that share the contig's node, since the
			// provider merges against that node. A contig is the coarsest safe flush boundary; the
			// buffer is also flushed mid-contig once it fills up (see handleStore()).
			if (batch != null && !batch.isEmpty()) {
				radixStore.updateBatch(batch, provider);
			}
		}

		@Override
		public boolean isAllowMoreKmers() {
			return true;
		}

		@Override
		protected boolean handleStore(long kmer) {
			if (batch != null) {
				if (batch.add(kmer)) {
					radixStore.updateBatch(batch, provider);
				}
				// The counted-k-mer return value is unused by this reader (endContig does no contig
				// bookkeeping and isAllowMoreKmers() is always true), so the deferred move result of a
				// batched k-mer need not be reported here.
				return false;
			}
			return store.update(kmer, provider);
		}
	}
}