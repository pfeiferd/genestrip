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
package org.metagene.genestrip.fastq;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.Thread.UncaughtExceptionHandler;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.io.BufferedLineReader;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.util.ByteArrayUtil;

public abstract class AbstractFastqReader {
	private static final byte[] LINE_3 = new byte[] { '+', '\n' };

	protected final Log logger = LogFactory.getLog("fastqreader");

	protected long reads;
	protected long kMers;

	private final BufferedLineReader bufferedLineReaderFastQ;
	private final ReadEntry[] readStructPool;
	private final BlockingQueue<ReadEntry> blockingQueue;
	private final Thread[] consumers;
	protected final int k;

	private boolean dump;
//	private Thread mainThread;
//	private boolean readsDone;

	private final List<Throwable> throwablesInThreads;

	public AbstractFastqReader(int k, int maxReadSizeBytes, int maxQueueSize, int consumerNumber, Object... config) {
		this.k = k;
		bufferedLineReaderFastQ = new BufferedLineReader();

		readStructPool = new ReadEntry[consumerNumber == 0 ? 1 : (maxQueueSize + consumerNumber + 1)];
		for (int i = 0; i < readStructPool.length; i++) {
			readStructPool[i] = createReadEntry(maxReadSizeBytes, config);
		}
		blockingQueue = consumerNumber == 0 ? null
				: new ArrayBlockingQueue<AbstractFastqReader.ReadEntry>(maxQueueSize);

		consumers = new Thread[consumerNumber];
		throwablesInThreads = Collections.synchronizedList(new ArrayList<Throwable>());

		for (int i = 0; i < consumers.length; i++) {
			consumers[i] = createAndStartThread(i, config);
		}
	}

	protected Thread createAndStartThread(int i, Object... config) {
		Thread t = new Thread(createRunnable(i, config));
		t.setName("Fastq reader thread #" + i);
		t.setUncaughtExceptionHandler(new UncaughtExceptionHandler() {
			@Override
			public void uncaughtException(Thread t, Throwable e) {
				throwablesInThreads.add(e);
			}
		});
		t.start();

		return t;
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

	protected Runnable createRunnable(int rindex, Object... config) {
		return new Runnable() {
			private int index = rindex;

			@Override
			public void run() {
				while (!dump) {
					try {
						ReadEntry readStruct = blockingQueue.take();
//						try {
							nextEntry(readStruct, index);
//						} finally {
//							if (readsDone) {
//								synchronized (mainThread) {
//									readStruct.pooled = true;
//									mainThread.notify();
//								}
//							} else {
								readStruct.pooled = true;
//							}
//						}
//					}
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

	protected ReadEntry createReadEntry(int maxReadSizeBytes, Object... config) {
		return new ReadEntry(maxReadSizeBytes);
	}

	public void dump() {
		dump = true;
		for (int i = 0; i < consumers.length; i++) {
			consumers[i].interrupt();
		}
	}

	protected Log getLogger() {
		return logger;
	}

	protected void readFastq(File file) throws IOException {
		if (logger.isInfoEnabled()) {
			logger.info("Reading fastq file " + file);
		}
		InputStream inputStream = StreamProvider.getInputStreamForFile(file);
		readFastq(inputStream);
		inputStream.close();
		if (logger.isInfoEnabled()) {
			logger.info("Closed fastq file " + file);
		}
	}

	protected void readFastq(InputStream inputStream) throws IOException {
		reads = 0;
		kMers = 0;
		bufferedLineReaderFastQ.setInputStream(inputStream);
//		readsDone = false;
//		mainThread = Thread.currentThread();

		start();

		if (logger.isInfoEnabled()) {
			logger.info("Number of consumer threads: " + consumers.length);
		}
		ReadEntry readStruct = nextFreeReadStruct();
		for (readStruct.readDescriptorSize = bufferedLineReaderFastQ.nextLine(readStruct.readDescriptor)
				- 1; readStruct.readDescriptorSize >= 0; readStruct.readDescriptorSize = bufferedLineReaderFastQ
						.nextLine(readStruct.readDescriptor) - 1) {
			readStruct.readDescriptor[readStruct.readDescriptorSize] = 0;
			readStruct.readSize = bufferedLineReaderFastQ.nextLine(readStruct.read) - 1;
			readStruct.read[readStruct.readSize] = 0;


// In some fastq files there occur two lines with '+' in a row which messes up the whole read process.
// This loop compensates for this bug in fastq files:
//			
//			readStruct.readProbsSize = bufferedLineReaderFastQ.nextLine(readStruct.readProbs) - 1;
//			readStruct.readProbs[readStruct.readProbsSize] = 0;
//			while (readStruct.readProbsSize == 1 && readStruct.readProbs[0] == '+') {
//				readStruct.readProbsSize = bufferedLineReaderFastQ.nextLine(readStruct.readProbs) - 1;
//				readStruct.readProbs[readStruct.readProbsSize] = 0;
//			}

			// Original code:			
			bufferedLineReaderFastQ.skipLine(); // Ignoring line 3.
			readStruct.readProbsSize = bufferedLineReaderFastQ.nextLine(readStruct.readProbs) - 1;
			readStruct.readProbs[readStruct.readProbsSize] = 0;
			
			readStruct.readNo = reads;

			reads++;
			kMers += readStruct.readSize - k + 1;
			if (blockingQueue == null) {
				nextEntry(readStruct, 0);
			} else {
				try {
					blockingQueue.put(readStruct);
				} catch (InterruptedException e) {
					throw new RuntimeException(e);
				}
			}
			checkAndLogConsumerThreadProblem();
			readStruct = nextFreeReadStruct();
			log();
		}
// This newer approach caused the threads to hang at the end of the reading process.
// I can't be bothered to find out why, so I use the older (ugly) polling approach which works well...		
//		readsDone = true;
//		if (blockingQueue != null) {
//			readStruct.pooled = true;
//			boolean stillWorking = true;
//			while (stillWorking) {
//				checkAndLogConsumerThreadProblem();
//				stillWorking = false;
//				for (int i = 0; i < readStructPool.length; i++) {
//					synchronized (mainThread) {
//						if (!readStructPool[i].pooled) {
//							stillWorking = true;
//							try {
//								mainThread.wait();
//							} catch (InterruptedException e) {
//								// Ignore.
//							}
//							break;
//						}
//					}
//				}
//			}
//		}
		if (blockingQueue != null) {
			readStruct.pooled = true;
			// Gentle polling and waiting until all consumers are done. 
			boolean stillWorking = true;
			while (stillWorking) {
				checkAndLogConsumerThreadProblem();
				stillWorking = false;
				for (int i = 0; i < readStructPool.length; i++) {
					if (!readStructPool[i].pooled) {
						stillWorking = true;
						try {
							Thread.sleep(100);
						} catch (InterruptedException e) {
							// Ignore.
						}
						break;
					}
				}
			}				
		}
		
		
		if (logger.isInfoEnabled()) {
			logger.info("Total number of reads: " + reads);
		}
		done();
	}

	protected void log() {
	}

	private ReadEntry nextFreeReadStruct() {
		if (blockingQueue == null) {
			return readStructPool[0];
		}
		for (int i = 0; i < readStructPool.length; i++) {
			if (readStructPool[i].pooled) {
				readStructPool[i].pooled = false;
				return readStructPool[i];
			}
		}
		throw new IllegalStateException("There should always be a read struct available...");
	}

	protected long getMillis() {
		return bufferedLineReaderFastQ.getMillis();
	}

	protected void rewriteInput(ReadEntry readStruct, OutputStream out) throws IOException {
		synchronized (out) {
			out.write(readStruct.readDescriptor, 0, readStruct.readDescriptorSize);
			out.write('\n');
			out.write(readStruct.read, 0, readStruct.readSize);
			out.write('\n');
			out.write(LINE_3);
			out.write(readStruct.readProbs, 0, readStruct.readProbsSize);
			out.write('\n');
			updateWriteStats();
		}
	}

	protected void updateWriteStats() {
	}

	// Must be thread safe. Can freely operate on readStruct.
	protected abstract void nextEntry(ReadEntry readStruct, int threadIndex) throws IOException;

	protected void done() throws IOException {
	};

	protected void start() throws IOException {
	};

	protected static class ReadEntry {
		public long readNo;
		public boolean pooled;

		public final byte[] readDescriptor;
		public int readDescriptorSize;

		public final byte[] read;
		public int readSize;

		public final byte[] readProbs;
		public int readProbsSize;

		protected ReadEntry(int maxReadSizeBytes) {
			readDescriptor = new byte[maxReadSizeBytes];
			read = new byte[maxReadSizeBytes];
			readProbs = new byte[maxReadSizeBytes];

			pooled = true;
		}

		public void print(PrintStream printStream) {
			ByteArrayUtil.println(readDescriptor, printStream);
			ByteArrayUtil.println(read, printStream);
			printStream.println('+');
			ByteArrayUtil.println(readProbs, printStream);
		}
	}
}
