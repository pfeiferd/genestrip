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
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
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
	private final ObjectGoal<AccessionMap, GSProject> accessionTrieGoal;
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final ObjectGoal<KMerStoreWrapper, GSProject> filledStoreGoal;

	private UpdatedStoreGoal updatedStoreGoal;

	private boolean dump;
	private final Thread[] consumers;
	private boolean[] working;
	private Object[] syncs = new Object[256];

	@SafeVarargs
	public UpdateStoreGoal(GSProject project, String name, ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal,
			ObjectGoal<TaxTree, GSProject> taxTreeGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
			ObjectGoal<AccessionMap, GSProject> accessionTrieGoal,
			ObjectGoal<KMerStoreWrapper, GSProject> filledStoreGoal, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.SER),
				Goal.append(deps, categoriesGoal, taxTreeGoal, fnaFilesGoal, accessionTrieGoal, filledStoreGoal));
		this.categoriesGoal = categoriesGoal;
		this.taxTreeGoal = taxTreeGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.accessionTrieGoal = accessionTrieGoal;
		this.filledStoreGoal = filledStoreGoal;
		consumers = new Thread[project.getConfig().getThreads()];
		working = new boolean[consumers.length];
	}

	public void setUpdatedStoreGoal(UpdatedStoreGoal updatedStoreGoal) {
		this.updatedStoreGoal = updatedStoreGoal;
	}

	@Override
	public void makeFile(File storeFile) {
		try {
			KMerStoreWrapper wrapper = filledStoreGoal.get();
			KMerSortedArray<String> store = wrapper.getKmerStore();
			AbstractFastaReader fastaReader = null;
			
			for (int i = 0; i < syncs.length; i++) {
				syncs[i] = new Object();
			}

			BlockingQueue<File> blockingQueue = new ArrayBlockingQueue<File>(consumers.length);
			for (int i = 0; i < consumers.length; i++) {
				consumers[i] = createAndStartThread(createFastaReaderRunnable(i, store, blockingQueue), i);
			}
			if (consumers.length == 0) {
				fastaReader = createFastaReader(store);
			}

			for (File fnaFile : fnaFilesGoal.getFiles()) {
				RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
				if (categoriesGoal.get()[0].contains(cat)) {
					if (fastaReader == null) {
						try {
							blockingQueue.put(fnaFile);
						} catch (InterruptedException e) {
							throw new RuntimeException(e);
						}
					} else {
						fastaReader.readFasta(fnaFile);
					}
				}
			}

			if (fastaReader == null) {
				// Gentle polling and waiting until all consumers are done.
				boolean stillWorking = true;
				while (stillWorking) {
					stillWorking = false;
					for (int i = 0; i < working.length; i++) {
						if (working[i]) {
							stillWorking = true;
							break;
						}
					}
					try {
						Thread.sleep(100);
					} catch (InterruptedException e) {
						// Ignore.
					}
				}
			}

			KMerStoreWrapper wrapper2 = new KMerStoreWrapper((KMerSortedArray<String>) store);
			wrapper2.save(storeFile);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File saved " + storeFile);
			}
			updatedStoreGoal.setStoreWrapper(wrapper2);
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			dump();
		}
	}

	protected Thread createAndStartThread(Runnable runnable, int i) {
		Thread t = new Thread(runnable);
		t.setName("Store updater thread #" + i);
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
						File fnaFile = blockingQueue.take();
						working[i] = true;
						fastaReader.readFasta(fnaFile);
						working[i] = false;
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

	protected AbstractFastaReader createFastaReader(KMerSortedArray<String> store) {
		return new MyFastaReader(getProject().getConfig().getMaxReadSizeBytes(), taxTreeGoal.get(),
				accessionTrieGoal.get(), store);
	}

	public void dump() {
		dump = true;
		for (int i = 0; i < consumers.length; i++) {
			consumers[i].interrupt();
		}
	}

	protected  class MyFastaReader extends AbstractStoreFastaReader {
		private final KMerSortedArray<String> store;
		private final UpdateValueProvider<String> provider;

		public MyFastaReader(int bufferSize, TaxTree taxTree, AccessionMap accessionMap,
				KMerSortedArray<String> store) {
			super(bufferSize, null, accessionMap, store.getK());
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
					return syncs[(int)(position % syncs.length)];
				}
			};
			includeRegion = true;
		}

		@Override
		protected void infoLine() throws IOException {
			byteRingBuffer.reset();
			updateNodeFromInfoLine();
		}

		@Override
		protected void handleStore() {
			// TODO: Synchronize me efficiently...
			store.update(byteRingBuffer, provider, false);
		}
	}
}