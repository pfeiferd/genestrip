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
import java.util.*;

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
import org.metagene.genestrip.store.KMerStore;
import org.metagene.genestrip.store.RadixKMerStore;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;

/**
 * Goal ({@code FILL_DB}) that allocates the k-mer store (a {@link RadixKMerStore} or a sorted array
 * sized from the estimated {@link FillBloomFilterGoal.DBSize}), fills it by labelling each k-mer
 * with its tax node's taxid, marks the requested nodes, optimizes the store and produces the
 * indexed {@link Database}.
 *
 * @param <P> the project type
 */
public class FillDBGoal<P extends GSProject> extends FastaReaderGoal<Database, P>  implements Goal.LogHeapInfo {
	private final ObjectGoal<AccessionMap, P> accessionMapGoal;
	private final ObjectGoal<FillBloomFilterGoal.DBSize, P> sizeGoal;
	private final ObjectGoal<TaxTree, P> taxTreeGoal;
	private final List<MyFastaReader> readers;

	private KMerStore<String> store;

	/**
	 * Creates the goal, wiring the tax tree, accession-map and estimated-size goals it uses to build
	 * and fill the store.
	 *
	 * @param project the project type
	 * @param bundle the execution context
	 * @param categoriesGoal the goal providing the RefSeq categories to include
	 * @param taxNodesGoal the goal providing the requested tax nodes
	 * @param taxTreeGoal the goal providing the taxonomy tree
	 * @param fnaFilesGoal the goal providing the downloaded fna files
	 * @param additionalGoal the goal providing additional files mapped to tax nodes
	 * @param accessionMapGoal the goal providing the accession map
	 * @param sizeGoal the goal providing the estimated database size
	 * @param deps additional goal dependencies
	 */
	@SafeVarargs
	public FillDBGoal(P project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
					  ObjectGoal<Set<TaxIdNode>, P> taxNodesGoal,
					  ObjectGoal<TaxTree, P> taxTreeGoal,
					  RefSeqFnaFilesDownloadGoal<P> fnaFilesGoal,
					  ObjectGoal<Map<File, TaxIdNode>, P> additionalGoal,
					  ObjectGoal<AccessionMap, P> accessionMapGoal,
					  ObjectGoal<FillBloomFilterGoal.DBSize, P> sizeGoal,
					  Goal<P>... deps) {
		super(project, GSGoalKey.FILL_DB, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, Goal.append(deps, taxTreeGoal, accessionMapGoal, sizeGoal));
		this.accessionMapGoal = accessionMapGoal;
		this.sizeGoal = sizeGoal;
		this.taxTreeGoal = taxTreeGoal;
		readers = new ArrayList<>();
	}

	@Override
	protected void doMakeThis() {
		int k = intConfigValue(GSConfigKey.KMER_SIZE);
		double fillFpp = doubleConfigValue(GSConfigKey.FILL_BLOOM_FILTER_FPP);
		double optFpp = doubleConfigValue(GSConfigKey.OPT_BLOOM_FILTER_FPP);
		boolean xor = booleanConfigValue(GSConfigKey.XOR_BLOOM_HASH);
		double resizeFactor = doubleConfigValue(GSConfigKey.DB_RESIZING_FACTOR);
		boolean useRadixStore = booleanConfigValue(GSConfigKey.USE_RADIX_STORE);
		FillBloomFilterGoal.DBSize dbSize = sizeGoal.get();

		if (useRadixStore) {
			// The radix store reserves capacity per radix bucket; the resizing factor scales each
			// bucket (analogous to scaling the total for the sorted array). The radix width must
			// match the one FillBloomFilterGoal used to count the per-bucket sizes.
			int radixBits = intConfigValue(GSConfigKey.RADIX_STORE_BITS);
			int[] bucketSizes = scaleBucketSizes(dbSize.getBucketSizes(), resizeFactor);
			store = new RadixKMerStore<>(k, radixBits, bucketSizes, fillFpp, optFpp, null, xor);
		} else {
			long dedupSize = dbSize.getSize();
			long size = resizeFactor == 1d ? dedupSize : (long) (dedupSize * resizeFactor);
			KMerSortedArray<String> array = new KMerSortedArray<>(k, fillFpp, optFpp, null, false, xor);
			array.initSize(size);
			store = array;
		}
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Store size in kmers: " + store.getSize());
			getLogger().info("DB Size in MB (without Bloom filter): " + (store.getSize() * 10) / (1024 * 1024));
		}

		try {
			readFastas();
			TaxTree taxTree = taxTreeGoal.get();
			// New nodes might have been added as file or id nodes. We must therefore reestablish
			// pre-order positions:
			taxTree.reinitPositions();
			long tooManyCounter = 0;
			for (MyFastaReader reader : readers) {
				tooManyCounter += reader.tooManyCounter;
			}
			if (getLogger().isWarnEnabled() && tooManyCounter > 0) {
				getLogger().warn("Not stored kmers: " + tooManyCounter);
			}
			long unused = store.getSize() - store.getEntries();
			if (getLogger().isInfoEnabled() && unused > 0) {
				getLogger().info("Unused kmers spots: " + unused +
						" (corresponds to " + ((100d * unused) / store.getSize()) + " %)");
			}
			SmallTaxTree smallTaxTree = taxTree.toSmallTaxTree();
			for (TaxIdNode node : taxNodesGoal.get()) {
				SmallTaxIdNode smallNode = smallTaxTree.getNodeByTaxId(node.getTaxId());
				if (smallNode != null) {
					smallNode.setRequested(true);
				}
			}
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Sorting kmers ...");
			}
			store.optimize();
			Database database = new Database(store, smallTaxTree, getProject().getAllAsProperties());
			database.initStoreIndices();
			set(database);
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			store = null;
			cleanUpThreads();
		}
	}

	// Scales the per-radix-bucket capacities by the resizing factor (a new array is returned; the
	// original from the size goal is not modified). resizeFactor == 1 returns the array unchanged.
	private static int[] scaleBucketSizes(int[] bucketSizes, double resizeFactor) {
		if (bucketSizes == null) {
			throw new IllegalStateException(
					"useRadixStore is set but the size goal did not provide per-bucket k-mer counts.");
		}
		if (resizeFactor == 1d) {
			return bucketSizes;
		}
		int[] scaled = new int[bucketSizes.length];
		for (int i = 0; i < bucketSizes.length; i++) {
			scaled[i] = (int) (bucketSizes[i] * resizeFactor);
		}
		return scaled;
	}

	@Override
	protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
		TaxTree taxTree = taxTreeGoal.get();
		boolean idNodes = booleanConfigValue(GSConfigKey.ID_NODES);
		boolean fileNodes = booleanConfigValue(GSConfigKey.FILE_NODES);
		boolean dataNodes = booleanConfigValue(GSConfigKey.DATA_NODES);

		byte[] idBuffer = new byte[128];
		idBuffer[0] = '0';
		idBuffer[1] = '0';

		TaxTree.IDStringGenerator idStringGenerator = new TaxTree.IDStringGenerator() {
			@Override
			public String generateID(int counter) {
				int len = ByteArrayUtil.intToByteArray(counter, idBuffer, 2);
				return new String(idBuffer, 0, len);
			}
		};

		MyFastaReader fastaReader = new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
				taxNodesGoal.get(), isIncludeRefSeqFna() ? accessionMapGoal.get() : null, store,
				intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
				(Rank) configValue(GSConfigKey.MAX_GENOMES_PER_TAXID_RANK),
				longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID),
				intConfigValue(GSConfigKey.MAX_DUST),
				intConfigValue(GSConfigKey.STEP_SIZE),
				booleanConfigValue(GSConfigKey.COMPLETE_GENOMES_ONLY),
				regionsPerTaxid,
				booleanConfigValue(GSConfigKey.ENABLE_LOWERCASE_BASES)) {
			@Override
			protected TaxIdNode reworkNode() {
				node.markRequired();
				TaxIdNode res = node;
				if (dataNodes) {
					if (Rank.DATA.ordinal() != res.getRankOrdinal()) {
						res = taxTree.dataNode(res, idStringGenerator);
					}
				}
				if (fileNodes && file != null) {
					if (Rank.FILE.ordinal() != res.getRankOrdinal()) {
						res = taxTree.fileNode(res, file.getName(), idStringGenerator);
					}
				}
				if (idNodes) {
					if (Rank.ID.ordinal() != res.getRankOrdinal()) {
						int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
						if (pos < 0) {
							// No space in the header: the id is the rest of the line. Trim the
							// trailing line terminator(s) so the id-node name excludes '\n' (and '\r').
							pos = size;
							while (pos > 0 && (target[pos - 1] == '\n' || target[pos - 1] == '\r')) {
								pos--;
							}
						}
						res = taxTree.idNode(res, target, 1, pos, idStringGenerator);
					}
				}
				res.markRequired();
				return res;
			}
		};
		readers.add(fastaReader);
		return fastaReader;
	}

	/**
	 * FASTA reader that stores each k-mer under its tax node's taxid, counting the k-mers dropped
	 * once the store is full.
	 */
	protected static class MyFastaReader extends AbstractStoreFastaReader {
		private final KMerStore<String> store;
		private long tooManyCounter;

		/**
		 * Creates the reader storing k-mers into the given store.
		 *
		 * @param bufferSize the read buffer size
		 * @param taxNodes the requested tax nodes
		 * @param accessionMap the accession-to-taxid map
		 * @param store the k-mer store to fill
		 * @param maxGenomesPerTaxId the maximum number of genomes per tax id
		 * @param maxGenomesPerTaxIdRank the rank at which the per-tax-id genome limit applies
		 * @param maxKmersPerTaxId the maximum number of k-mers per tax id
		 * @param maxDust the maximum allowed low-complexity (dust) run length
		 * @param stepSize the k-mer sampling step size
		 * @param completeGenomesOnly whether to include only complete genomes
		 * @param regionsPerTaxid the per-taxid region trie
		 * @param enableLowerCaseBases whether lower-case bases are included
		 */
		public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap,
							 KMerStore<String> store, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases) {
			super(bufferSize, taxNodes, accessionMap, store.getK(), maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, maxDust, stepSize, completeGenomesOnly, regionsPerTaxid, enableLowerCaseBases);
			this.store = store;
		}

		@Override
		protected boolean handleStore() {
			if (store.isFull()) {
				tooManyCounter++;
			} else {
				if (node.getTaxId() != null) {
					return store.putLong(byteRingBuffer.getStandardKMer(), node.getTaxId());
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