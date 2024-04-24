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
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import org.apache.commons.logging.Log;
import org.metagene.genestrip.io.BufferedLineReader;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.ExecutorServiceBundle;
import org.metagene.genestrip.util.GSLogFactory;

public abstract class AbstractFastqReader {
	private static final byte[] LINE_3 = new byte[] { '+', '\n' };

	protected final Log logger = GSLogFactory.getLog("fastqreader");

	protected long reads;
	protected long kMers;

	private final BufferedLineReader bufferedLineReaderFastQ;
	private final ReadEntry[] readStructPool;
	private final BlockingQueue<ReadEntry> blockingQueue;
	protected final int k;
	private final byte[] plusLine;

	private boolean dump;

	protected final ExecutorServiceBundle bundle;

	public AbstractFastqReader(int k, int initialSizeBytes, int maxQueueSize, ExecutorServiceBundle bundle,
			Object... config) {
		this.k = k;
		this.bundle = bundle;
		this.plusLine = new byte[initialSizeBytes];
		bufferedLineReaderFastQ = new BufferedLineReader();
		int consumerNumber = bundle.getThreads();

		readStructPool = new ReadEntry[consumerNumber == 0 ? 1 : (maxQueueSize + consumerNumber + 1)];
		for (int i = 0; i < readStructPool.length; i++) {
			readStructPool[i] = createReadEntry(initialSizeBytes, config);
		}
		blockingQueue = consumerNumber == 0 ? null
				: new ArrayBlockingQueue<AbstractFastqReader.ReadEntry>(maxQueueSize);

		if (blockingQueue != null) {
			for (int i = 0; i < consumerNumber; i++) {
				bundle.execute(createRunnable(i, config));
			}
		}
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

	protected ReadEntry createReadEntry(int initialReadSizeBytes, Object... config) {
		return new ReadEntry(initialReadSizeBytes);
	}

	public void dump() {
		dump = true;
		bundle.interruptAll();
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

	protected void growReadBuffer(ReadEntry readStruct, BufferedLineReader lineReader) throws IOException {
		while (readStruct.readSize == readStruct.read.length) {
			byte[] newBuffer = new byte[readStruct.read.length * 2];
			System.arraycopy(readStruct.read, 0, newBuffer, 0, readStruct.read.length);
			readStruct.readSize = lineReader.nextLine(newBuffer, readStruct.read.length) - 1;
			readStruct.read = newBuffer;
		}
		// Probs must grow in the same way.
		readStruct.readProbs = new byte[readStruct.read.length];
	}

	protected void readFastq(InputStream inputStream) throws IOException {
		reads = 0;
		kMers = 0;
		bufferedLineReaderFastQ.setInputStream(inputStream);

		start();

		if (logger.isInfoEnabled()) {
			logger.info("Number of consumer threads: " + bundle.getThreads());
		}
		ReadEntry readStruct = nextFreeReadStruct();
		for (readStruct.readDescriptorSize = bufferedLineReaderFastQ.nextLine(readStruct.readDescriptor)
				- 1; readStruct.readDescriptorSize >= 0; readStruct.readDescriptorSize = bufferedLineReaderFastQ
						.nextLine(readStruct.readDescriptor) - 1) {
			readStruct.readDescriptor[readStruct.readDescriptorSize] = 0;
			readStruct.readSize = bufferedLineReaderFastQ.nextLine(readStruct.read) - 1;
			if (readStruct.readSize == readStruct.read.length) {
				growReadBuffer(readStruct, bufferedLineReaderFastQ);
			}
			readStruct.read[readStruct.readSize] = 0;

			// Old code:
			// bufferedLineReaderFastQ.skipLine(); // Ignoring line 3.
			// New code does basic consistency check:
			bufferedLineReaderFastQ.nextLine(plusLine);
			if (plusLine[0] != '+') {
				throw new IllegalStateException("Inconsistent read no. " + reads + ", the line with '+' is missing.");
			}

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
					throw new FastqReaderInterruptedException(e);
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
	}

	protected void start() throws IOException {
	}

	protected static class ReadEntry {
		public long readNo;
		public boolean pooled;

		public final byte[] readDescriptor;
		public int readDescriptorSize;

		public byte[] read;
		public int readSize;

		public byte[] readProbs;
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
