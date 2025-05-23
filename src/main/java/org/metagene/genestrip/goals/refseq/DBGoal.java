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
import java.time.temporal.ChronoUnit;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import me.tongfei.progressbar.DelegatingProgressBarConsumer;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarBuilder;
import me.tongfei.progressbar.ProgressBarStyle;
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
import org.metagene.genestrip.util.GSLogFactory;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;

public class DBGoal extends ObjectGoal<Database, GSProject> {
	private final ObjectGoal<Set<RefSeqCategory>, GSProject> categoriesGoal;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionTrieGoal;
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final ObjectGoal<Database, GSProject> filledStoreGoal;
	private final ExecutionContext bundle;

	private boolean dump;
	private int doneCounter;
	private Object[] syncs = new Object[256];

	private boolean minUpdate;
	private ProgressBar progressBar;

	@SafeVarargs
	public DBGoal(GSProject project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, GSProject> categoriesGoal,
				  ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
				  ObjectGoal<TaxTree, GSProject> taxTreeGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
				  ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalGoal,
			ObjectGoal<AccessionMap, GSProject> accessionTrieGoal, ObjectGoal<Database, GSProject> filledStoreGoal,
			Goal<GSProject>... deps) {
		super(project, GSGoalKey.UPDATE_DB, Goal.append(deps, additionalGoal, categoriesGoal, taxTreeGoal, taxNodesGoal,
				fnaFilesGoal, accessionTrieGoal, filledStoreGoal));
		this.categoriesGoal = categoriesGoal;
		this.taxNodesGoal = taxNodesGoal;
		this.taxTreeGoal = taxTreeGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.additionalGoal = additionalGoal;
		this.accessionTrieGoal = accessionTrieGoal;
		this.filledStoreGoal = filledStoreGoal;
		this.bundle = bundle;
		minUpdate = project.booleanConfigValue(GSConfigKey.MIN_UPDATE);
	}
	
	@Override
	protected void doMakeThis() {       
		try {
			Database wrapper = filledStoreGoal.get();
			KMerSortedArray<String> store = wrapper.getKmerStore();
			AbstractRefSeqFastaReader fastaReader = createFastaReader(store);

			for (int i = 0; i < syncs.length; i++) {
				syncs[i] = new Object();
			}
			bundle.clearThrowableList();

			BlockingQueue<FileAndNode> blockingQueue = null;
			// For minUpdate the fna files must be read in the same order as during the goals from before.
			// Otherwise, the wrong k-mers will be compared to the ones from the DB.
			// For !minUpdate multi-threading can be enabled.
			if (bundle.getThreads() > 0) {
				blockingQueue = new ArrayBlockingQueue<FileAndNode>(intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE));
				for (int i = 0; i < bundle.getThreads(); i++) {
					bundle.execute(createFastaReaderRunnable(i, store, blockingQueue));
				}
			}

			int sumFiles = 0;
			List<File> refSeqFiles = fnaFilesGoal.getFiles();
			sumFiles += refSeqFiles.size();
			Map<File, TaxTree.TaxIdNode> additionalMap = additionalGoal.get();
			sumFiles += additionalMap.size();
			try (ProgressBar pb = (progressBar = createProgressBar(sumFiles))) {
				doneCounter = 0;
				for (File fnaFile : refSeqFiles) {
					RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
					if (categoriesGoal.get().contains(cat)) {
						if (blockingQueue == null) {
							fastaReader.readFasta(fnaFile);
						} else {
							try {
								doneCounter++;
								blockingQueue.put(new FileAndNode(fnaFile, null));
							} catch (InterruptedException e) {
								throw new RuntimeException(e);
							}
						}
					}
					checkAndLogConsumerThreadProblem();
				}
				for (File additionalFasta : additionalMap.keySet()) {
					if (blockingQueue == null) {
						fastaReader.ignoreAccessionMap(additionalMap.get(additionalFasta));
						fastaReader.readFasta(additionalFasta);
					} else {
						try {
							doneCounter++;
							blockingQueue.put(new FileAndNode(additionalFasta, additionalMap.get(additionalFasta)));
						} catch (InterruptedException e) {
							throw new RuntimeException(e);
						}
					}
					checkAndLogConsumerThreadProblem();
				}
				// Gentle polling and waiting until all consumers are done.
				while (doneCounter > 0) {
					checkAndLogConsumerThreadProblem();
					try {
						Thread.sleep(100);
					} catch (InterruptedException e) {
						// Ignore.
					}
				}
			}
			store.fix();
			set(new Database(store, wrapper.getTaxTree()));
			if (getLogger().isTraceEnabled()) {
				getLogger().trace("KMers moved: " + store.getKMersMoved());
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			cleanUpThreads();
		}
	}

	protected ProgressBar createProgressBar(int max) {
		return booleanConfigValue(GSConfigKey.PROGRESS_BAR) ?
				GSProgressBarCreator.newGSProgressBar(getKey().getName(), max, 60000, " files", null, getLogger()) :
				null;
	}

	protected void checkAndLogConsumerThreadProblem() {
		if (!bundle.getThrowableList().isEmpty()) {
			for (Throwable t : bundle.getThrowableList()) {
				if (getLogger().isErrorEnabled()) {
					getLogger().error("Error in consumer thread: ", t);
				}
			}
			bundle.clearThrowableList();
			throw new RuntimeException("Error(s) in consumer thread(s).");
		}
	}

	protected Runnable createFastaReaderRunnable(int i, KMerSortedArray<String> store,
			BlockingQueue<FileAndNode> blockingQueue) {
		return new Runnable() {
			@Override
			public void run() {
				AbstractRefSeqFastaReader fastaReader = createFastaReader(store);

				while (!dump) {
					try {
						try {
							FileAndNode fileAndNode = blockingQueue.take();
							fastaReader.ignoreAccessionMap(fileAndNode.getNode());
							fastaReader.readFasta(fileAndNode.getFile());
							if (progressBar != null) {
								progressBar.step();
							}
						} finally {
							doneCounter--;
						}
					} catch (IOException e) {
						throw new RuntimeException(e);
					} catch (InterruptedException e) {
						if (!dump) {
							throw new RuntimeException(e);
						}
					}
				}
			}
		};
	}

	protected AbstractRefSeqFastaReader createFastaReader(KMerSortedArray<String> store) {
		return new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES), taxTreeGoal.get(), taxNodesGoal.get(),
				accessionTrieGoal.get(), store, intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
				(Rank) configValue(GSConfigKey.MAX_GENOMES_PER_TAXID_RANK),
				longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID),
				intConfigValue(GSConfigKey.MAX_DUST),
				intConfigValue(GSConfigKey.STEP_SIZE),
				booleanConfigValue(GSConfigKey.UPDATE_WITH_COMPLETE_GENOMES_ONLY));
	}

	public void dump() {
		super.dump();
		cleanUpThreads();
	}
	
	protected void cleanUpThreads() {
		dump = true;
		bundle.interruptAll();		
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

				@Override
				public Object getSynchronizationObject(long position) {
					return syncs[(int) (position % syncs.length)];
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
			return store.update(byteRingBuffer.getKMer(), byteRingBuffer.getReverseKMer(), provider);
		}
	}

	protected static final class FileAndNode {
		private final File file;
		private final TaxIdNode node;

		public FileAndNode(File file, TaxIdNode node) {
			this.file = file;
			this.node = node;
		}

		public File getFile() {
			return file;
		}

		public TaxIdNode getNode() {
			return node;
		}
	}
}