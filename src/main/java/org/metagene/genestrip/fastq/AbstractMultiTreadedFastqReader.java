package org.metagene.genestrip.fastq;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.util.BufferedLineReader;
import org.metagene.genestrip.util.StreamProvider;

public abstract class AbstractMultiTreadedFastqReader {
	// private static final byte[] LINE_3 = new byte[] { '+', '\n' };

	protected final Log logger = LogFactory.getLog(getClass());

	protected long reads;
	private final BufferedLineReader bufferedLineReaderFastQ;
	private final ReadEntry[] readStructPool;
	private final BlockingQueue<ReadEntry> blockingQueue;
	private final Thread[] consumers;
	
	private boolean closed;
	
	public AbstractMultiTreadedFastqReader(int maxReadSizeBytes, int maxQueueSize, int consumerNumber) {
		bufferedLineReaderFastQ = new BufferedLineReader();

		readStructPool = new ReadEntry[maxQueueSize + consumerNumber + 1];
		for (int i = 0; i < readStructPool.length; i++) {
			readStructPool[i] = new ReadEntry(maxReadSizeBytes);
		}
		blockingQueue = new ArrayBlockingQueue<AbstractMultiTreadedFastqReader.ReadEntry>(maxQueueSize);

		Runnable runnable = new Runnable() {
			@Override
			public void run() {
				while (!closed) {
					try {
						ReadEntry readStruct = blockingQueue.take();
						nextEntry(readStruct);
						readStruct.pooled = true;
					} catch (IOException e) {
						throw new RuntimeException(e);
					} catch (InterruptedException e) {
						if (!closed) {
							throw new RuntimeException(e);
						}
					}
				}
			}
		};
		
		consumers = new Thread[consumerNumber];
		for (int i = 0; i < consumers.length; i++) {
			consumers[i] = new Thread(runnable);
			consumers[i].start();
		}
	}
	
	public void close() {
		closed = true;
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
		bufferedLineReaderFastQ.setInputStream(inputStream);

		start();

		ReadEntry readStruct = nextFreeReadStruct();
		for (readStruct.readDescriptorSize = bufferedLineReaderFastQ.nextLine(readStruct.readDescriptor)
				- 1; readStruct.readDescriptorSize > 0; readStruct.readDescriptorSize = bufferedLineReaderFastQ
						.nextLine(readStruct.readDescriptor) - 1) {
			readStruct.readDescriptor[readStruct.readDescriptorSize] = 0;
			readStruct.readSize = bufferedLineReaderFastQ.nextLine(readStruct.read) - 1;
			readStruct.read[readStruct.readSize] = 0;
			bufferedLineReaderFastQ.skipLine(); // Ignoring line 3.
			readStruct.readProbsSize = bufferedLineReaderFastQ.nextLine(readStruct.readProbs) - 1;
			readStruct.readProbs[readStruct.readProbsSize] = 0;

			reads++;
			try {
				blockingQueue.put(readStruct);
			} catch (InterruptedException e) {
				throw new RuntimeException(e);
			}
			readStruct = nextFreeReadStruct();
		}
		if (logger.isInfoEnabled()) {
			logger.info("Total number of reads: " + reads);
		}
		done();
	}

	private ReadEntry nextFreeReadStruct() {
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

	/*
	 * protected void rewriteInput(OutputStream out) throws IOException {
	 * out.write(readDescriptor, 0, readDescriptorSize); out.write('\n');
	 * out.write(read, 0, readSize); out.write('\n'); out.write(LINE_3);
	 * out.write(readProbs, 0, readProbsSize); out.write('\n'); }
	 */

	// Must be thread safe. Can freely operate on readStruct.
	protected abstract void nextEntry(ReadEntry readStruct) throws IOException;

	protected void done() throws IOException {
	};

	protected void start() throws IOException {
	};

	public static class ReadEntry {
		public boolean pooled;

		public final byte[] readDescriptor;
		public int readDescriptorSize;

		public final byte[] read;
		public int readSize;

		public final byte[] readProbs;
		public int readProbsSize;

		private ReadEntry(int maxReadSizeBytes) {
			readDescriptor = new byte[maxReadSizeBytes];
			read = new byte[maxReadSizeBytes];
			readProbs = new byte[maxReadSizeBytes];

			pooled = true;
		}
	}
}
