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
import java.lang.Thread.UncaughtExceptionHandler;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerSortedArray.UpdateValueProvider;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class UpdateStoreGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionTrieGoal;
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final ObjectGoal<KMerStoreWrapper, GSProject> filledStoreGoal;

	private UpdatedStoreGoal updatedStoreGoal;

	private boolean dump;
	private final Thread[] consumers;
	private int doneCounter;
	private Object[] syncs = new Object[256];

	private final List<Throwable> throwablesInThreads;

	@SafeVarargs
	public UpdateStoreGoal(GSProject project, String name, ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal,
			ObjectGoal<TaxTree, GSProject> taxTreeGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
			ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalGoal,
			ObjectGoal<AccessionMap, GSProject> accessionTrieGoal,
			ObjectGoal<KMerStoreWrapper, GSProject> filledStoreGoal, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.DB),
				Goal.append(deps, categoriesGoal, taxTreeGoal, fnaFilesGoal, accessionTrieGoal, filledStoreGoal));
		this.categoriesGoal = categoriesGoal;
		this.taxTreeGoal = taxTreeGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.additionalGoal = additionalGoal;
		this.accessionTrieGoal = accessionTrieGoal;
		this.filledStoreGoal = filledStoreGoal;
		consumers = new Thread[project.getConfig().getThreads()];
		throwablesInThreads = Collections.synchronizedList(new ArrayList<Throwable>());
	}

	public void setUpdatedStoreGoal(UpdatedStoreGoal updatedStoreGoal) {
		this.updatedStoreGoal = updatedStoreGoal;
	}

	@Override
	public void makeFile(File storeFile) {
		try {
			KMerStoreWrapper wrapper = filledStoreGoal.get();
			KMerSortedArray<String> store = wrapper.getKmerStore();
			AbstractRefSeqFastaReader fastaReader = createFastaReader(store);

			for (int i = 0; i < syncs.length; i++) {
				syncs[i] = new Object();
			}
			throwablesInThreads.clear();

			BlockingQueue<File> blockingQueue = new ArrayBlockingQueue<File>(consumers.length);
			for (int i = 0; i < consumers.length; i++) {
				consumers[i] = createAndStartThread(createFastaReaderRunnable(i, store, blockingQueue), i);
			}

			doneCounter = 0;
			for (File fnaFile : fnaFilesGoal.getFiles()) {
				RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
				if (categoriesGoal.get()[1].contains(cat)) {
					if (consumers.length == 0) {
						fastaReader.readFasta(fnaFile);
					} else {
						try {
							doneCounter++;
							blockingQueue.put(fnaFile);
						} catch (InterruptedException e) {
							throw new RuntimeException(e);
						}
					}
				}
				checkAndLogConsumerThreadProblem();
			}
			// Simply do this on main thread - should not be so much work...
			Map<File, TaxIdNode> additionalMap = additionalGoal.get();
			for (File additionalFasta : additionalMap.keySet()) {
				fastaReader.ignoreAccessionMap(additionalMap.get(additionalFasta));
				fastaReader.readFasta(additionalFasta);
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

			KMerStoreWrapper wrapper2 = new KMerStoreWrapper((KMerSortedArray<String>) store);
			File filterFile = getProject().getFilterFile(this);
			wrapper2.save(storeFile, filterFile);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File saved " + storeFile + " along with index " + filterFile);
			}
			updatedStoreGoal.setStoreWrapper(wrapper2);
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			dump();
		}
	}

	protected void checkAndLogConsumerThreadProblem() {
		if (!throwablesInThreads.isEmpty()) {
			for (Throwable t : throwablesInThreads) {
				if (getLogger().isErrorEnabled()) {
					getLogger().error("Error in consumer thread: ", t);
				}
			}
			throw new RuntimeException("Error(s) in consumer thread(s).");
		}
	}

	protected Thread createAndStartThread(Runnable runnable, int i) {
		Thread t = new Thread(runnable);
		t.setName("Store updater thread #" + i);
		t.setUncaughtExceptionHandler(new UncaughtExceptionHandler() {
			@Override
			public void uncaughtException(Thread t, Throwable e) {
				throwablesInThreads.add(e);
			}
		});
		t.start();

		return t;
	}

	protected Runnable createFastaReaderRunnable(int i, KMerSortedArray<String> store,
			BlockingQueue<File> blockingQueue) {
		return new Runnable() {
			@Override
			public void run() {
				AbstractFastaReader fastaReader = createFastaReader(store);

				while (!dump) {
					try {
						try {
							File fnaFile = blockingQueue.take();
							fastaReader.readFasta(fnaFile);
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
		return new MyFastaReader(getProject().getConfig().getMaxReadSizeBytes(), taxTreeGoal.get(),
				accessionTrieGoal.get(), store, getProject().getMaxGenomesPerTaxid(), getProject().getMaxDust());
	}

	public void dump() {
		dump = true;
		for (int i = 0; i < consumers.length; i++) {
			consumers[i].interrupt();
		}
	}

	protected class MyFastaReader extends AbstractStoreFastaReader {
		private final KMerSortedArray<String> store;
		private final UpdateValueProvider<String> provider;

		public MyFastaReader(int bufferSize, TaxTree taxTree, AccessionMap accessionMap, KMerSortedArray<String> store,
				int maxGenomesPerTaxId, int maxDust) {
			super(bufferSize, null, accessionMap, store.getK(), maxGenomesPerTaxId, maxDust);
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
			includeRegion = true;
		}

		@Override
		protected void infoLine() throws IOException {
			byteRingBuffer.reset();
			if (!ignoreMap) {
				updateNodeFromInfoLine();
			}
		}

		@Override
		protected void handleStore() {
			store.update(byteRingBuffer.getKMer(), provider);
		}
	}
}