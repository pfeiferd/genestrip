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
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;

public class FillDBGoal<P extends GSProject> extends FastaReaderGoal<Database, P>  implements Goal.LogHeapInfo {
	private final ObjectGoal<AccessionMap, P> accessionMapGoal;
	private final ObjectGoal<Long, P> sizeGoal;
	private final ObjectGoal<TaxTree, P> taxTreeGoal;
	private final List<MyFastaReader> readers;

	private KMerSortedArray<String> store;

	@SafeVarargs
	public FillDBGoal(P project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
					  ObjectGoal<Set<TaxIdNode>, P> taxNodesGoal,
					  ObjectGoal<TaxTree, P> taxTreeGoal,
					  RefSeqFnaFilesDownloadGoal fnaFilesGoal,
					  ObjectGoal<Map<File, TaxIdNode>, P> additionalGoal,
					  ObjectGoal<AccessionMap, P> accessionMapGoal,
					  ObjectGoal<Long, P> sizeGoal,
					  Goal<P>... deps) {
		super(project, GSGoalKey.FILL_DB, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, Goal.append(deps, taxTreeGoal, accessionMapGoal, sizeGoal));
		this.accessionMapGoal = accessionMapGoal;
		this.sizeGoal = sizeGoal;
		this.taxTreeGoal = taxTreeGoal;
		readers = new ArrayList<>();
	}

	@Override
	protected void doMakeThis() {
		store = new KMerSortedArray<>(intConfigValue(GSConfigKey.KMER_SIZE),
				doubleConfigValue(GSConfigKey.BLOOM_FILTER_FPP), doubleConfigValue(GSConfigKey.OPT_BLOOM_FILTER_FPP), null, false, booleanConfigValue(GSConfigKey.XOR_BLOOM_HASH));
		double resizeFactor = doubleConfigValue(GSConfigKey.DB_RESIZING_FACTOR);
		long size = resizeFactor == 1d ? sizeGoal.get() : (long) (sizeGoal.get() * resizeFactor);
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Store size in kmers: " + size);
			getLogger().info("DB Size in MB (without Bloom filter): " + (size * 10) / (1024 * 1024));
		}
		store.initSize(size);

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
			ensureAllTreeNodesInDB(smallTaxTree.getRoot(), store);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Sorting kmers ...");
			}
			store.optimize();
			Database wrapper = new Database(store, smallTaxTree);
			set(wrapper);
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			store = null;
			cleanUpThreads();
		}
	}

	@Override
	protected AbstractStoreFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
		TaxTree taxTree = taxTreeGoal.get();
		boolean idNodes = booleanConfigValue(GSConfigKey.ID_NODES);
		boolean fileNodes = booleanConfigValue(GSConfigKey.FILE_NODES);

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
				if (fileNodes && file != null) {
					if (!Rank.FILE.equals(res.getRank())) {
						res = taxTree.fileNode(res, file.getName(), idStringGenerator);
					}
				}
				if (idNodes) {
					if (!Rank.ID.equals(res.getRank())) {
						int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
						if (pos < 0) {
							pos = size;
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
							 KMerSortedArray<String> store, int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int maxDust, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean enableLowerCaseBases) {
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