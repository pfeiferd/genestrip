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
package org.metagene.genestrip.match;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.fastq.AbstractLoggingFastqStreamer;
import org.metagene.genestrip.io.*;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerUniqueCounter;
import org.metagene.genestrip.store.KMerUniqueCounterBits;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CGAT;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class FastqKMerMatcher extends AbstractLoggingFastqStreamer {
    private final static SmallTaxIdNode INVALID_NODE = new SmallTaxIdNode("INVALID");

    protected final KMerSortedArray<SmallTaxIdNode> kmerStore;
    protected final int maxKmerResCounts;

    // Turned from KMerUniqueCounter to KMerUniqueCounterBits for potential method inlining.
    protected KMerUniqueCounterBits uniqueCounter;
    protected final CountsPerTaxid[] statsIndex;
    protected final long readNoPerCPerStat[][];

    protected final int maxPaths;
    protected final SmallTaxTree taxTree;
    protected final double maxReadTaxErrorCount;
    protected final double maxReadClassErrorCount;
    protected OutputStream indexed;

    // This should stay a box type for the line root.get(taxid.getTaxId(),
    // maxReadSize);
    protected final Integer initialReadSize;

    // A PrintStream is implicitly synchronized. So we don't need to worry about
    // multi-threading when using it.
    protected PrintStream out;
    protected final boolean writeAll;

    protected int consumers;

    private AfterMatchCallback afterMatchCallback;

    public FastqKMerMatcher(KMerSortedArray<SmallTaxIdNode> kmerStore, int initialReadSize, int maxQueueSize,
                            ExecutionContext bundle, boolean withProbs, int maxKmerResCounts, SmallTaxTree taxTree, int maxPaths,
                            double maxReadTaxErrorCount, double maxReadClassErrorCount, boolean writeAll) {
        super(kmerStore.getK(), initialReadSize, maxQueueSize, bundle, withProbs, maxPaths);
        consumers = bundle.getThreads() <= 0 ? 1 : bundle.getThreads();
        this.kmerStore = kmerStore;
        this.statsIndex = new CountsPerTaxid[kmerStore.getNValues()];
        this.readNoPerCPerStat = new long[consumers][kmerStore.getNValues()];
        this.initialReadSize = initialReadSize;
        this.maxKmerResCounts = maxKmerResCounts;
        this.taxTree = taxTree;
        this.maxReadTaxErrorCount = maxReadTaxErrorCount;
        this.maxReadClassErrorCount = maxReadClassErrorCount;
        this.maxPaths = maxPaths;
        this.writeAll = writeAll;
        if (taxTree != null) {
            taxTree.initCountSize(consumers);
        }
    }

    @Override
    protected ReadEntry createReadEntry(int initialReadSizeBytes, boolean withProbs, Object... config) {
        return new MyReadEntry(initialReadSizeBytes, withProbs, (int) config[0]);
    }

    public MatchingResult runMatcher(StreamingResource fastq, File filteredFile, File krakenOutStyleFile,
                                     KMerUniqueCounter uniqueCounter) throws IOException {
        return runMatcher(new StreamingResourceListStream(fastq), filteredFile, krakenOutStyleFile, uniqueCounter);
    }

    public MatchingResult runMatcher(StreamingResourceStream fastqs, File filteredFile, File krakenOutStyleFile,
                                     KMerUniqueCounter uniqueCounter) throws IOException {
        try (OutputStream lindexed = filteredFile != null ? StreamProvider.getOutputStreamForFile(filteredFile) : null;
             // A PrintStream is implicitly synchronized. So we don't need to worry about
             // multi threading when using it.
             PrintStream lout = krakenOutStyleFile != null
                     ? new PrintStream(StreamProvider.getOutputStreamForFile(krakenOutStyleFile))
                     : null) {
            indexed = lindexed;
            out = lout;

            initStats();
            initUniqueCounter(uniqueCounter);
            processFastqStreams(fastqs);
        }
        out = null;
        indexed = null;

        Map<String, CountsPerTaxid> taxid2Stats = new HashMap<>();
        for (CountsPerTaxid stats : statsIndex) {
            if (stats != null) {
                taxid2Stats.put(stats.getTaxid(), stats);
            }
        }

        Map<String, short[]> countMap = null;
        if (uniqueCounter != null) {
            Object2LongMap<String> counts = uniqueCounter.getUniqueKmerCounts();
            for (CountsPerTaxid stats : statsIndex) {
                if (stats != null) {
                    stats.uniqueKmers = counts.getLong(stats.getTaxid());
                }
            }
            if (uniqueCounter instanceof KMerUniqueCounterBits) {
                if (((KMerUniqueCounterBits) uniqueCounter).isWithCounts()) {
                    countMap = ((KMerUniqueCounterBits) uniqueCounter).getMaxCountsCounts(maxKmerResCounts);
                    for (CountsPerTaxid stats : statsIndex) {
                        if (stats != null) {
                            stats.maxKMerCounts = countMap.get(stats.getTaxid());
                        }
                    }
                }
            }
            this.uniqueCounter = null;
        } else {
            for (CountsPerTaxid stats : statsIndex) {
                if (stats != null) {
                    stats.uniqueKmers = -1;
                }
            }
        }

        return new MatchingResult(kmerStore.getK(), taxid2Stats, totalReads, totalKMers, totalBPs,
                countMap == null ? null : countMap.get(null));
    }

    @Override
    protected void readFastq(InputStream inputStream) throws IOException {
        try {
            if (taxTree != null) {
                taxTree.resetCounts(this);
            }
            for (long[] a : readNoPerCPerStat) {
                Arrays.fill(a, -1);
            }
            super.readFastq(inputStream);
        } finally {
            if (taxTree != null) {
                taxTree.releaseOwner();
            }
        }
    }

    // Package private for testing purposes.
    void initStats() {
        Arrays.fill(statsIndex, null);
    }

    void initUniqueCounter(KMerUniqueCounter uniqueCounter) {
        // Turned from KMerUniqueCounter to KMerUniqueCounterBits for potential method inlining.
        this.uniqueCounter = (KMerUniqueCounterBits) uniqueCounter;
        if (uniqueCounter != null) {
            uniqueCounter.clear();
        }
    }

    @Override
    // Made final for potential inlining by JVM
    protected final void nextEntry(ReadEntry entry, int index) throws IOException {
        MyReadEntry myEntry = (MyReadEntry) entry;
        myEntry.bufferPos = 0;

        myEntry.usedPaths = 0;
        myEntry.classNode = null;
        for (int i = 0; i < maxPaths; i++) {
            myEntry.readTaxIdNode[i] = null;
            myEntry.counts[i] = 0;
        }

        boolean found =  matchRead(myEntry, index);
        afterMatch(myEntry, found);
        if (afterMatchCallback != null) {
            afterMatchCallback.afterMatch(myEntry, found);
        }
    }

    public void setAfterMatchCallback(AfterMatchCallback afterMatchCallback) {
        this.afterMatchCallback = afterMatchCallback;
    }

    protected void afterMatch(MyReadEntry myEntry, boolean found) throws IOException {
        if (found && indexed != null) {
            rewriteInput(myEntry, indexed);
        }
        if (out != null) {
            if (writeAll || myEntry.classNode != null) {
                synchronized (out) {
                    myEntry.flush(out);
                }
            }
        }
    }

    private byte[] kmerHelp = new byte[31];

    // Made final for potential inlining by JVM
    protected final boolean matchRead(final MyReadEntry entry, final int index) {
        boolean found = false;
        int prints = 0;
        int readTaxErrorCount = taxTree == null ? -1 : 0;

        SmallTaxIdNode taxIdNode;
        int max = entry.readSize - k + 1;
        SmallTaxIdNode lastTaxid = null;
        int contigLen = 0;
        CountsPerTaxid stats = null;

        long kmer = -1;
        long reverseKmer = -1;
        int oldIndex = 0;
        for (int i = 0; i < max; i++) {
            if (kmer == -1) {
                kmer = CGAT.kMerToLongStraight(entry.read, i, k, entry.badPos);
                if (kmer == -1) {
                    oldIndex = i;
                    i = entry.badPos[0];
                } else {
                    reverseKmer = CGAT.kMerToLongReverse(entry.read, i, k, null);
                }
            } else {
                kmer = CGAT.nextKMerStraight(kmer, entry.read[i + k - 1], k);
                if (kmer == -1) {
                    oldIndex = i;
                    i += k - 1;
                } else {
                    reverseKmer = CGAT.nextKMerReverse(reverseKmer, entry.read[i + k - 1], k);
                }
            }
            taxIdNode = kmer == -1 ? INVALID_NODE :
                    kmerStore.getLong(CGAT.standardKMer(kmer, reverseKmer), entry.indexPos);
            if (readTaxErrorCount != -1) {
                if (taxIdNode == null || taxIdNode == INVALID_NODE) {
                    readTaxErrorCount++;
                    if (maxReadTaxErrorCount >= 0) {
                        if ((maxReadTaxErrorCount >= 1 && readTaxErrorCount > maxReadTaxErrorCount)
                                || (readTaxErrorCount > maxReadTaxErrorCount * max)) {
                            readTaxErrorCount = -1;
                        }
                    }
                } else {
                    updateReadTaxid(taxIdNode, entry, index);
                }
            }
            if (taxIdNode != lastTaxid) {
                if (contigLen > 0) {
                    if (out != null) {
                        printKrakenStyleOut(entry, lastTaxid, contigLen, prints++);
                    }
                    if (stats != null) {
                        synchronized (stats) {
                            stats.contigs++;
                            stats.contigLenSquaredSum += contigLen * contigLen;
                            if (contigLen > stats.maxContigLen) {
                                stats.maxContigLen = contigLen;
                                int j = 1;
                                for (; j < entry.readDescriptorSize && entry.readDescriptor[j] != ' '; j++) {
                                    stats.maxContigDescriptor[j - 1] = entry.readDescriptor[j];
                                }
                                stats.maxContigDescriptor[j - 1] = 0;
                            }
                        }
                    }
                    contigLen = 0;
                }
            }
            if (taxIdNode == INVALID_NODE) {
                contigLen += i >= max ? max - oldIndex : i - oldIndex + 1;
            }
            else {
                contigLen++;
            }
            lastTaxid = taxIdNode;
            if (taxIdNode != null && taxIdNode != INVALID_NODE) {
                short vi = taxIdNode.getStoreIndex();
                stats = getCountsPerTaxid(taxIdNode, vi);
                synchronized (stats) {
                    stats.kmers++;
                    found = true;
                    if (readNoPerCPerStat[index][vi] != entry.readNo) {
                        stats.reads1KMer++;
                        readNoPerCPerStat[index][vi] = entry.readNo;
                    }
                }
                if (uniqueCounter != null) {
                    uniqueCounter.put(CGAT.standardKMer(kmer, reverseKmer), taxIdNode.getTaxId(), entry.indexPos[0]);
                }
            } else {
                stats = null;
            }
        }
        if (contigLen > 0 && out != null) {
            printKrakenStyleOut(entry, lastTaxid, contigLen, prints);
        }
        if (found) {
            if (contigLen > 0) {
                if (stats != null) {
                    synchronized (stats) {
                        stats.contigs++;
                        stats.contigLenSquaredSum += contigLen * contigLen;
                        if (contigLen > stats.maxContigLen) {
                            stats.maxContigLen = contigLen;
                            int j = 1;
                            for (; j < entry.readDescriptorSize && entry.readDescriptor[j] != ' '; j++) {
                                stats.maxContigDescriptor[j - 1] = entry.readDescriptor[j];
                            }
                            stats.maxContigDescriptor[j - 1] = 0;
                        }
                    }
                }
            }
            if (readTaxErrorCount != -1) {
                int ties = 0;
                for (int i = 0; i < entry.usedPaths; i++) {
                    short sum = taxTree.sumCounts(entry.readTaxIdNode[i], index, entry.readNo);
                    if (sum > entry.counts[0]) {
                        entry.counts[0] = sum;
                        entry.readTaxIdNode[0] = entry.readTaxIdNode[i];
                        ties = 0;
                    } else if (sum == entry.counts[0]) {
                        ties++;
                        entry.counts[ties] = sum;
                        entry.readTaxIdNode[ties] = entry.readTaxIdNode[i];
                    }
                }
                SmallTaxIdNode node = entry.readTaxIdNode[0];
                for (int i = 1; i <= ties; i++) {
                    node = taxTree.getLeastCommonAncestor(node, entry.readTaxIdNode[i]);
                }
                entry.classNode = node;
                // For 'readKmers', I decided to count in the k-mers from 'entry.readTaxIdNode[0]' and not just 'node'.
                // (They only differ in case of a tie anyways.) But if there is tie, then the k-mers from one of the tie's nodes
                // solidify the LCA in a sense - so the counts from one the involved paths are included.
                int readKmers = ties > 0 ? taxTree.sumCounts(entry.readTaxIdNode[0], index, entry.readNo) : entry.counts[0];
                int classErrC = max - readKmers;
                if (maxReadClassErrorCount < 0 || (maxReadClassErrorCount >= 1 && classErrC <= maxReadClassErrorCount)
                        || (classErrC <= maxReadClassErrorCount * max)) {
                    double err = ((double) readTaxErrorCount) / max;
                    double classErr = ((double) classErrC) / max;
                    entry.classNode = node;
                    short vi = node.getStoreIndex();
                    if (vi != -1) {
                        stats = getCountsPerTaxid(node, vi);
                        synchronized (stats) {
                            stats.reads++;
                            stats.readsKmers += readKmers;
                            stats.readsBPs += entry.readSize;
                            stats.errorSum += err;
                            stats.errorSquaredSum += err * err;
                            stats.classErrorSum += classErr;
                            stats.classErrorSquaredSum += classErr * classErr;
                        }
                    }
                    else if (getLogger().isWarnEnabled()) {
                        getLogger().warn("Missing database entry for tax node: " + node);
                    }
                }
            }
        }

        return found;
    }

    protected final CountsPerTaxid getCountsPerTaxid(final SmallTaxIdNode node, final short vi) {
        CountsPerTaxid stats = statsIndex[vi];
        if (stats == null) {
            synchronized (statsIndex) {
                if (statsIndex[vi] == null) {
                    statsIndex[vi] = new CountsPerTaxid(node.getTaxId(), initialReadSize);
                }
                stats = statsIndex[vi];
            }
        }
        return stats;
    }

    // Made final for potential inlining by JVM
    protected final void updateReadTaxid(final SmallTaxIdNode node, final MyReadEntry entry, final int index) {
        taxTree.incCount(node, index, entry.readNo);

        boolean found = false;
        for (int i = 0; i < entry.usedPaths; i++) {
            if (taxTree.isAncestorOf(node, entry.readTaxIdNode[i])) {
                entry.readTaxIdNode[i] = node;
                found = true;
                break;
            } else if (taxTree.isAncestorOf(entry.readTaxIdNode[i], node)) {
                found = true;
                break;
            }
        }
        if (!found) {
            if (entry.usedPaths < maxPaths) {
                entry.readTaxIdNode[entry.usedPaths] = node;
                entry.usedPaths++;
            }
        }
    }

    protected void printKrakenStyleOut(final MyReadEntry entry, final SmallTaxIdNode taxid, final int contigLen, final int state) {
        if (state != 0) {
            entry.printChar(' ');
        }
        if (taxid == INVALID_NODE) {
            entry.printChar('A');
        }
        else if (taxid == null) {
            entry.printChar('0');
        } else {
            entry.printString(taxid.getTaxId());
        }
        entry.printChar(':');
        entry.printInt(contigLen);
    }

    public static class MyReadEntry extends ReadEntry {
        public byte[] buffer;
        public int bufferPos;
        public int[] badPos = new int[1];

        public int usedPaths;
        public SmallTaxIdNode[] readTaxIdNode;
        public short[] counts;
        public long[] indexPos;
        public SmallTaxIdNode classNode;

        public MyReadEntry(int maxReadSizeBytes, boolean withProbs, int paths) {
            super(maxReadSizeBytes, withProbs);

            buffer = null;
            readTaxIdNode = new SmallTaxIdNode[paths];
            counts = new short[paths];
            indexPos = new long[1];
        }

        public void printChar(final char c) {
            growPrintBuffer(1);
            buffer[bufferPos++] = (byte) c;
        }

        public void printString(final String s) {
            int len = s.length();
            growPrintBuffer(len);
            for (int i = 0; i < len; i++) {
                buffer[bufferPos++] = (byte) s.charAt(i);
            }
        }

        protected void growPrintBuffer(int additionalSize) {
            if (buffer == null) {
                buffer = new byte[additionalSize + 1024];
                return;
            }

            int newLen = buffer.length;
            while (bufferPos + additionalSize > newLen) {
                newLen *= 2;
            }
            if (newLen > buffer.length) {
                byte[] newBuffer = new byte[newLen];
                System.arraycopy(buffer, 0, newBuffer, 0, bufferPos);
                buffer = newBuffer;
            }
        }

        public void printInt(final int value) {
            growPrintBuffer(11); // Decimal version of int hat 11 bytes max.
            bufferPos = ByteArrayUtil.intToByteArray(value, buffer, bufferPos);
        }

        public void flush(PrintStream out) {
            if (buffer == null) {
                return;
            }
            if (classNode == null) {
                out.print("U\t");
            } else {
                out.print("C\t");
            }
            int index = ByteArrayUtil.indexOf(readDescriptor, 1, readDescriptor.length, ' ');
            ByteArrayUtil.print(readDescriptor, 1, index == -1 ? readDescriptor.length : index, out);
            out.print('\t');
            if (classNode == null) {
                out.print('0');
            } else {
                out.print(classNode.getTaxId());
            }
            out.print('\t');
            out.print(readSize);
            out.print('\t');
            out.write(buffer, 0, bufferPos);
            out.println();
        }
   }

   public interface AfterMatchCallback {
        void afterMatch(MyReadEntry entry, boolean found);
   }
}
