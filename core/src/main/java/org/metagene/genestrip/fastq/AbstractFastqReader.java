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

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;

import org.apache.commons.logging.Log;
import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.io.BufferedLineReader;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.GSLogFactory;
import org.metagene.genestrip.util.SimpleBlockingQueue;

/**
 * Abstract multi-threaded FASTQ/FASTA reader. It parses reads into a pool of reusable
 * {@link ReadEntry} buffers and either processes them inline or hands them to a configurable number
 * of consumer threads via a blocking queue, calling {@link #nextEntry(ReadEntry, int)} for each.
 */
public abstract class AbstractFastqReader {
    private static final byte[] LINE_3 = new byte[]{'+', '\n'};

    /** The logger used by this reader. */
    protected final Log logger = GSLogFactory.getLog("fastqreader");

    /** Number of reads processed so far. */
    protected long reads;
    /** Number of k-mers counted across all reads processed so far. */
    protected long kMers;
    /** Total number of base pairs across all reads processed so far. */
    protected long readBPs;

    private final BufferedLineReader bufferedLineReaderFastQ;
    private final ReadEntry[] readStructPool;
    private final BlockingQueue<ReadEntry> blockingQueue;
    private final byte[] plusLine;
    /** The k-mer length used when counting k-mers per read. */
    protected final int k;

    // volatile: written by dump() and read by the producer loop on a different thread.
    private volatile boolean dump;

    /** The execution context providing the consumer threads and error handling. */
    protected final ExecutionContext bundle;

    /**
     * Creates a reader for k-mers of length {@code k}, allocating read buffers of
     * {@code initialSizeBytes} and starting the number of consumer threads given by the execution
     * context (or processing inline if that number is zero).
     *
     * @param k the k-mer length used when counting k-mers per read.
     * @param initialSizeBytes the initial size in bytes of each allocated read buffer.
     * @param maxQueueSize the maximum number of read entries buffered in the blocking queue.
     * @param bundle the execution context providing the consumer threads and error handling.
     * @param withProbs whether reads carry per-base quality (probability) values.
     * @param config optional configuration passed on to the read-entry and runnable factories.
     */
    public AbstractFastqReader(int k, int initialSizeBytes, int maxQueueSize, ExecutionContext bundle,
                               boolean withProbs, Object... config) {
        this.k = k;
        this.bundle = bundle;
        bufferedLineReaderFastQ = new BufferedLineReader();
        int consumerNumber = bundle.getThreads();
        plusLine = new byte[initialSizeBytes];
        readStructPool = new ReadEntry[consumerNumber == 0 ? 2 : (maxQueueSize + consumerNumber + 2)];
        for (int i = 0; i < readStructPool.length; i++) {
            readStructPool[i] = createReadEntry(initialSizeBytes, withProbs, config);
        }
        blockingQueue = consumerNumber == 0 ? null
                : createBlockingQueue(maxQueueSize);

        if (blockingQueue != null) {
            for (int i = 0; i < consumerNumber; i++) {
                bundle.execute(createRunnable(i, config));
            }
        }
    }

    /**
     * Creates the blocking queue used to hand read entries to the consumer threads.
     *
     * @param maxQueueSize the maximum number of read entries the queue may hold.
     * @return the blocking queue used to pass read entries to the consumer threads.
     */
    protected BlockingQueue<ReadEntry> createBlockingQueue(int maxQueueSize) {
        // This simple blocking queue gives about 5% to 10% performance boost (on my Mac)
        // over the ArrayBlockingQueue. Also, it does not cause any memory churn
        // (unlike ArrayBlockingQueue).
        return new SimpleBlockingQueue<>(maxQueueSize);
        // return new ArrayBlockingQueue<>(maxQueueSize);
    }

    /**
     * Logs any exceptions raised on consumer threads and rethrows the first as a runtime exception.
     */
    // Made final for potential inlining by JVM
    protected final void checkAndLogConsumerThreadProblem() {
        synchronized (bundle) {
            if (bundle.hasThrowables()) {
                // To avoid co modification exception due to mo
                List<Throwable> b = new ArrayList<>(bundle.getThrowableList());
                for (Throwable t : b) {
                    if (getLogger().isErrorEnabled()) {
                        getLogger().error("Error in consumer thread: ", t);
                    }
                }
                bundle.clearThrowableList();
                if (b.size() > 0) {
                    // Let's take the first one of the potentially  many exceptions and wrap it.
                    throw new RuntimeException("Error(s) in consumer thread(s).", b.get(0));
                } else {
                    throw new RuntimeException("Error(s) in consumer thread(s).");
                }
            }
        }
    }

    /**
     * Creates the worker {@link Runnable} for a consumer thread, which takes read entries from the
     * queue, processes them via {@link #nextEntry(ReadEntry, int)} and returns them to the pool.
     *
     * @param rindex the index identifying this consumer thread.
     * @param config optional configuration for the consumer thread.
     * @return the runnable executed by the consumer thread.
     */
    protected Runnable createRunnable(int rindex, Object... config) {
        return new Runnable() {
            private final int index = rindex;

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

    /**
     * Creates a pooled read-entry buffer; subclasses may override to use a specialized entry.
     *
     * @param initialReadSizeBytes the initial size in bytes of the entry's buffers.
     * @param withProbs whether the entry carries per-base quality (probability) values.
     * @param config optional configuration for the created entry.
     * @return the newly created read entry.
     */
    protected ReadEntry createReadEntry(int initialReadSizeBytes, boolean withProbs, Object... config) {
        return new ReadEntry(initialReadSizeBytes, withProbs);
    }

    /**
     * Signals the reader to stop and interrupts all consumer threads.
     */
    public void dump() {
        dump = true;
        bundle.interruptAll();
    }

    /**
     * Returns the logger used by this reader.
     *
     * @return the logger used by this reader.
     */
    protected Log getLogger() {
        return logger;
    }

    /**
     * Reads all reads from the given stream (as FASTA if {@code fasta} is true, otherwise FASTQ),
     * dispatching them and waiting for all consumer threads to finish.
     *
     * @param inputStream the stream to read reads from.
     * @param fasta if {@code true} the input is parsed as FASTA, otherwise as FASTQ.
     * @throws IOException if the stream cannot be read.
     */
    protected void readFastq(InputStream inputStream, boolean fasta) throws IOException {
        start();
        reads = 0;
        kMers = 0;
        readBPs = 0;
        bufferedLineReaderFastQ.setInputStream(inputStream);

        if (fasta) {
            doReadFasta();
        } else {
            doReadFastq();
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

        done();
    }

    /**
     * Parses the input as FASTQ, dispatching one {@link ReadEntry} per read.
     *
     * @throws IOException if the input cannot be read.
     */
    protected void doReadFastq() throws IOException {
        ReadEntry readStruct = nextFreeReadStruct();
        for (readStruct.readDescriptorSize = bufferedLineReaderFastQ.nextLine(readStruct.readDescriptor)
                - 1; readStruct.readDescriptorSize >= 0; readStruct.readDescriptorSize = bufferedLineReaderFastQ
                .nextLine(readStruct.readDescriptor) - 1) {
            readStruct.readDescriptor[readStruct.readDescriptorSize] = 0;
            readStruct.readSize = bufferedLineReaderFastQ.nextLine(readStruct.read) - 1;
            if (readStruct.readSize == readStruct.read.length) {
                readStruct.growReadBuffer(bufferedLineReaderFastQ);
            }

            int newSize = bufferedLineReaderFastQ.nextLine(readStruct.read, readStruct.readSize) - 1;
            // If there is no '+' then the line still belongs to the read.
            while (readStruct.read[readStruct.readSize] != '+') {
                readStruct.readSize = newSize;
                if (newSize == readStruct.read.length) {
                    readStruct.growReadBuffer(bufferedLineReaderFastQ);
                }
                newSize = bufferedLineReaderFastQ.nextLine(readStruct.read, readStruct.readSize) - 1;
            }
            readStruct.read[readStruct.readSize] = 0;
            // This means the buffer (i.e. readStruct.read) was full but we found the '+',
            // but still not passed the plus line
            // as not read passed the buffer
            if (newSize == readStruct.read.length) {
                // Dump the rest of the plus line here (until eol reached).
                bufferedLineReaderFastQ.nextLine(plusLine);
            }

            // We got passed the '+' line...
            if (readStruct.readProbs != null) {
                int readProbsSize = bufferedLineReaderFastQ.nextLine(readStruct.readProbs) - 1;
                while (readProbsSize < readStruct.readSize) {
                    int oldSize = readProbsSize;
                    readProbsSize = bufferedLineReaderFastQ.nextLine(readStruct.readProbs, readProbsSize) - 1;
                    // The means we reached EOF:
                    if (readProbsSize == oldSize - 1) {
                        break;
                    }
                }
                readStruct.readProbsSize = readProbsSize;
                readStruct.readProbs[readProbsSize] = 0;
            } else {
                int readProbsSize = bufferedLineReaderFastQ.skipLine() - 1;
                while (readProbsSize < readStruct.readSize) {
                    int oldSize = readProbsSize;
                    readProbsSize = readProbsSize + bufferedLineReaderFastQ.skipLine() - 1;
                    // The means we reached EOF:
                    if (readProbsSize == oldSize - 1) {
                        break;
                    }
                }
                readStruct.readProbsSize = readProbsSize;
            }

            readStruct.readNo = reads;

            reads++;
            if (readStruct.readSize >= k) {
                kMers += readStruct.readSize - k + 1;
            }
            readBPs += readStruct.readSize;
            if (blockingQueue == null) {
                nextEntry(readStruct, 0);
                readStruct.pooled = true;
                if (dump) {
                    throw new FastqReaderInterruptedException();
                }
            } else {
                try {
                    blockingQueue.put(readStruct);
                } catch (InterruptedException e) {
                    throw new FastqReaderInterruptedException(e);
                }
            }
            checkAndLogConsumerThreadProblem();
            readStruct = nextFreeReadStruct();
            updateProgress();
        }
        readStruct.pooled = true;
    }

    /**
     * Parses the input as FASTA, dispatching one {@link ReadEntry} per sequence.
     *
     * @throws IOException if the input cannot be read.
     */
    protected void doReadFasta() throws IOException {
        ReadEntry readStruct = nextFreeReadStruct();
        for (readStruct.readDescriptorSize = bufferedLineReaderFastQ.nextLine(readStruct.readDescriptor)
                - 1; readStruct != null && readStruct.readDescriptorSize >= 0; ) {
            readStruct.readDescriptor[readStruct.readDescriptorSize] = 0;
            readStruct.readDescriptor[0] = '@';

            // ByteArrayUtil.println(readStruct.readDescriptor, System.out);

            readStruct.readSize = 0;
            int newSize = bufferedLineReaderFastQ.nextLine(readStruct.read, readStruct.readSize) - 1;
            // If there is no '>' then the line still belongs to the read.
            while (readStruct.read[readStruct.readSize] != '>' && newSize > readStruct.readSize - 1) {
                readStruct.readSize = newSize;
                if (newSize == readStruct.read.length) {
                    readStruct.growReadBuffer(bufferedLineReaderFastQ);
                }
                newSize = bufferedLineReaderFastQ.nextLine(readStruct.read, readStruct.readSize) - 1;
            }

            ReadEntry readStruct2 = null;
            // Not reached EOF?
            if (newSize != readStruct.readSize - 1) {
                readStruct2 = nextFreeReadStruct();
                int len = newSize - readStruct.readSize;
                System.arraycopy(readStruct.read, readStruct.readSize, readStruct2.readDescriptor, 0, len);
                if (newSize == readStruct.read.length) {
                    readStruct2.readDescriptorSize = bufferedLineReaderFastQ.nextLine(readStruct2.readDescriptor, len);
                } else {
                    readStruct2.readDescriptorSize = len;
                }
            }

            readStruct.read[readStruct.readSize] = 0;
//			ByteArrayUtil.println(readStruct.read, System.out);

            readStruct.readNo = reads;
            reads++;
            if (readStruct.readSize >= k) {
                kMers += readStruct.readSize - k + 1;
            }
            readBPs += readStruct.readSize;
            readStruct.readProbsSize = -1; // Indicate no probs available.
            if (blockingQueue == null) {
                nextEntry(readStruct, 0);
                readStruct.pooled = true;
                if (dump) {
                    throw new FastqReaderInterruptedException();
                }
            } else {
                try {
                    blockingQueue.put(readStruct);
                } catch (InterruptedException e) {
                    throw new FastqReaderInterruptedException(e);
                }
            }
            checkAndLogConsumerThreadProblem();
            updateProgress();
            readStruct = readStruct2;
        }
        if (readStruct != null) {
            readStruct.pooled = true;
        }
    }

    /**
     * Hook called after each read is dispatched so subclasses can report progress; does nothing by
     * default.
     */
    protected void updateProgress() {
    }

    private ReadEntry nextFreeReadStruct() {
        for (int i = 0; i < readStructPool.length; i++) {
            if (readStructPool[i].pooled) {
                readStructPool[i].pooled = false;
                return readStructPool[i];
            }
        }
        throw new IllegalStateException("There should always be a read structs available...");
    }

    /**
     * Writes the given read entry to the output stream in a thread-safe manner and updates the write
     * statistics.
     *
     * @param readStruct the read entry to write.
     * @param out the output stream to write to.
     * @throws IOException if writing to the stream fails.
     */
    protected void rewriteInput(ReadEntry readStruct, OutputStream out) throws IOException {
        synchronized (out) {
            readStruct.write(out);
            updateWriteStats();
        }
    }

    /**
     * Hook called after a read entry is written out so subclasses can update write statistics; does
     * nothing by default.
     */
    protected void updateWriteStats() {
    }

    /**
     * Processes a single read entry. Called on consumer threads, so it must be thread-safe; it may
     * freely operate on the given read entry.
     *
     * @param readStruct the read entry to process.
     * @param threadIndex the index of the consumer thread invoking this method.
     * @throws IOException if processing the read entry fails.
     */
    // Must be thread safe. Can freely operate on readStruct.
    protected abstract void nextEntry(ReadEntry readStruct, int threadIndex) throws IOException;

    /**
     * Hook called after all reads have been processed; does nothing by default.
     *
     * @throws IOException if the subclass fails to finalize processing.
     */
    protected void done() throws IOException {
    }

    /**
     * Hook called before reading begins; does nothing by default.
     *
     * @throws IOException if the subclass fails to initialize.
     */
    protected void start() throws IOException {
    }

    /**
     * A reusable, pooled buffer holding one read's descriptor, sequence and (optional) quality
     * bytes.
     */
    protected static class ReadEntry {
        /** The zero-based sequential number of this read. */
        public long readNo;
        /** Whether this entry is currently free in the pool and available for reuse. */
        // volatile: written by consumer threads, read by the producer's pool/poll loop.
        public volatile boolean pooled;

        /** The buffer holding the read's descriptor (header) bytes. */
        public final byte[] readDescriptor;
        /** The number of valid bytes in {@link #readDescriptor}. */
        public int readDescriptorSize;

        /** The buffer holding the read's sequence bytes. */
        public byte[] read;
        /** The number of valid bytes in {@link #read}. */
        public int readSize;

        /** The buffer holding the read's per-base quality bytes, or {@code null} if none. */
        public byte[] readProbs;
        /** The number of valid bytes in {@link #readProbs}, or {@code -1} if none. */
        public int readProbsSize;

        /**
         * Creates a read entry with buffers of the given size.
         *
         * @param maxReadSizeBytes the initial size in bytes of the entry's buffers.
         * @param withProbs whether a per-base quality buffer is allocated.
         */
        protected ReadEntry(int maxReadSizeBytes, boolean withProbs) {
            readDescriptor = new byte[maxReadSizeBytes];
            read = new byte[maxReadSizeBytes];
            readProbs = withProbs ? new byte[maxReadSizeBytes] : null;

            pooled = true;
        }

		/*
		public void print(PrintStream printStream) {
			ByteArrayUtil.println(readDescriptor, printStream);
			ByteArrayUtil.println(read, printStream);
			printStream.println('+');
			if (readProbs != null && readProbsSize >= 0) {
				ByteArrayUtil.println(readProbs, printStream);
			}
			else {
				for (int i = 0; i < readProbsSize; i++) {
					printStream.print('~');
				}
				printStream.println();
			}
		}
		 */

        /**
         * Writes this read to the stream in FASTQ format, synthesizing '~' quality values when none
         * are present.
         *
         * @param out the output stream to write to.
         * @throws IOException if writing to the stream fails.
         */
        public void write(OutputStream out) throws IOException {
            out.write(readDescriptor, 0, readDescriptorSize);
            out.write('\n');
            out.write(read, 0, readSize);
            out.write('\n');
            out.write(LINE_3);
            if (readProbs != null && readProbsSize >= 0) {
                out.write(readProbs, 0, readProbsSize);
            } else {
                for (int i = 0; i < readSize; i++) {
                    out.write('~');
                }
            }
            out.write('\n');
        }

        /**
         * Doubles the read buffer (and the probability buffer, if present) until the current line
         * fits, continuing to read the line from the given reader.
         *
         * @param lineReader the reader to continue reading the current line from.
         * @throws IOException if reading from the reader fails.
         */
        protected final void growReadBuffer(BufferedLineReader lineReader) throws IOException {
            while (readSize == read.length) {
                byte[] newBuffer = new byte[read.length * 2];
                System.arraycopy(read, 0, newBuffer, 0, read.length);
                readSize = lineReader.nextLine(newBuffer, read.length) - 1;
                read = newBuffer;
            }
            if (readProbs != null) {
                // Probs must grow in the same way.
                readProbs = new byte[read.length];
            }
        }
    }
}
