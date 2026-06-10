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

public abstract class AbstractFastqReader {
    private static final byte[] LINE_3 = new byte[]{'+', '\n'};

    protected final Log logger = GSLogFactory.getLog("fastqreader");

    protected long reads;
    protected long kMers;
    protected long readBPs;

    private final BufferedLineReader bufferedLineReaderFastQ;
    private final ReadEntry[] readStructPool;
    private final BlockingQueue<ReadEntry> blockingQueue;
    private final byte[] plusLine;
    protected final int k;

    private boolean dump;

    protected final ExecutionContext bundle;

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

    protected BlockingQueue<ReadEntry> createBlockingQueue(int maxQueueSize) {
        // This simple blocking queue gives about 5% to 10% performance boost (on my Mac)
        // over the ArrayBlockingQueue. Also, it does not cause any memory churn
        // (unlike ArrayBlockingQueue).
        return new SimpleBlockingQueue<>(maxQueueSize);
        // return new ArrayBlockingQueue<>(maxQueueSize);
    }

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

    protected ReadEntry createReadEntry(int initialReadSizeBytes, boolean withProbs, Object... config) {
        return new ReadEntry(initialReadSizeBytes, withProbs);
    }

    public void dump() {
        dump = true;
        bundle.interruptAll();
    }

    protected Log getLogger() {
        return logger;
    }

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
            if (readStruct.readSize > k) {
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
            if (readStruct.readSize > k) {
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

    protected void rewriteInput(ReadEntry readStruct, OutputStream out) throws IOException {
        synchronized (out) {
            readStruct.write(out);
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
