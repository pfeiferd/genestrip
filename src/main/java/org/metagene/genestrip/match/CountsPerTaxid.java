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

import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.Serializable;

public class CountsPerTaxid implements Serializable, Comparable<CountsPerTaxid> {
    public enum ValueType {
        READS("reads"), KMERS("kmers"), READS_BPS("reads bps"), READS_1KMER("read >=1 kmer"), READS_KMERS("reads kmers");

        private static final ValueType[] VALUES = ValueType.values();

        private final String name;

        ValueType(String name) {
            this.name = name;
        }

        public String getName() {
            return name;
        }
    }

    public static class AccValues implements Serializable {
        private long accumulated;
        private double accumulatedNormalized;

        public AccValues(long value, long dbKMers) {
            this.accumulated = value;
            this.accumulatedNormalized = dbKMers > 0 ? ((double) value) / dbKMers : 0;
        }

        public double getAccumulatedNormalized() {
            return accumulatedNormalized;
        }

        public long getAccumulated() {
            return accumulated;
        }

        public void accumulateFrom(AccValues extendedValues) {
            accumulated += extendedValues.getAccumulated();
            accumulatedNormalized += extendedValues.getAccumulatedNormalized();
        }
    }

    private static final long serialVersionUID = 1L;

    protected String taxid;
    protected long reads;
    protected long reads1KMer;
    protected long readsBPs;
    protected long readsKmers;
    protected long uniqueKmers; // All unique kmers counted - even from unclassified reads.
    protected long kmers;
    protected int contigs;
    protected int maxContigLen;
    protected byte[] maxContigDescriptor;
    protected short[] maxKMerCounts;
    protected double errorSum;
    protected double errorSquaredSum;
    protected double accErrorSum;
    protected double accErrorSquaredSum;

    // To complete values
    private int pos;
    private String name;
    private Rank rank;
    private long dbKMers;
    private String parentTaxId;
    private final AccValues[] extendedValues;

    public CountsPerTaxid(String taxid, int maxReadSizeBytes) {
        this.taxid = taxid;
        maxContigDescriptor = new byte[maxReadSizeBytes];
        extendedValues = new AccValues[ValueType.VALUES.length];
    }

    public CountsPerTaxid(String taxid, long totalReads, long totalKmers, short[] totalMaxCounts) {
        this.taxid = taxid;
        this.reads = totalReads;
        this.kmers = totalKmers;
        this.maxKMerCounts = totalMaxCounts;
        maxContigDescriptor = new byte[0];
        extendedValues = new AccValues[ValueType.VALUES.length];
    }

    @Override
    public int compareTo(CountsPerTaxid o) {
        return this.pos - o.pos;
    }

    @MDCDescription(pos = 0, name = "pos", desc = "Sort position of entry.")
    public int getPos() {
        return pos;
    }

    @MDCDescription(pos = 1, name = "name", desc = "The name associated with the tax id.")
    public String getName() {
        return name;
    }

    @MDCDescription(pos = 2, name= "rank", desc = "The rank of the tax id.")
    public Rank getRank() {
        return rank;
    }

    @MDCDescription(pos = 3, name = "taxid", desc = "The tax id.")
    public String getTaxid() {
        return taxid;
    }

    // Only classified if a read is ENTIRELY consistent with the kmers of a taxon (and/or its ancestors) in the DB
    @MDCDescription(pos = 4, name = "reads", desc = "The number of reads classified with respect to the tax id.")
    public long getReads() {
        return reads;
    }

    @MDCDescription(pos = 5, name="kmers from reads", desc ="The number of *k*-mers from classified reads which are consistent with the read's tax id.")
    public long getReadsKMers() {
        return readsKmers;
    }

    @MDCDescription(pos = 6, name="kmers", desc = "*All* matched *k*-mers which are specific to the tax id's genome (according to the database). The *k*-mers do not have to be in an accordingly classified read for this.")
    public long getKMers() {
        return kmers;
    }

    @MDCDescription(pos = 7, name="unique kmers", desc = "*All* unique *k*-mers, which are specific to the tax id's genome (according to the database). " +
            "Here, multiple occurrences of the same *k*-mer are only counted once. " +
            "Genestrip always performs exact counting according to [KrakenUniq's exact counting](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1568-0#Sec8). (Genestrip implements an efficient in-memory storage method for related counts based on a bit vectors.)")
    public long getUniqueKMers() {
        return uniqueKmers;
    }

    @MDCDescription(pos = 8, name="contigs", desc = "The number of contiguous sequences of *k*-mers that are specific to the tax id's genome.")
    public int getContigs() {
        return contigs;
    }

    @MDCDescription(pos = 9, name="average contig length", desc = "The average length of contiguous sequences of *k*-mers (as the number of *k*-mers) that are specific to the tax id's genome.")
    public double getAverageContigLen() {
        return ((double) kmers) / contigs;
    }

    @MDCDescription(pos = 10, name="max contig length", desc ="The maximum length of all contiguous sequences of *k*-mers (as the number of *k*-mers) that are specific to the tax id's genome.")
    public int getMaxContigLen() {
        return maxContigLen;
    }

    @MDCDescription(pos = 11, name = "reads >=1 kmer", desc = "Reads with at least on *k*-mer of the respective tax id.")
    public long getReads1KMer() {
        return reads1KMer;
    }

    @MDCDescription(pos = 12, name = "reads bps", desc = "The total number of base pairs of reads classified with respect to the tax id.")
    public long getReadsBPs() {
        return readsBPs;
    }

    @MDCDescription(pos = 13, name = "avg. read length", desc = "The average length of classified reads in base pairs.")
    public double getAverageReadLength() {
        return ((double) readsBPs) / reads;
    }

    @MDCDescription(pos = 14, name = "db coverage", desc ="The ratio `unique kmers` / u<sub>t</sub>, , where *u<sub>t</sub>* = `db kmers`")
    public double getCoverage() {
        return ((double) uniqueKmers) / dbKMers;
    }

    @MDCDescription(pos = 15, name = "exp. unique kmers", desc ="The number of expected unique *k*-mers, which is *u<sub>t</sub> * (1 - (1 - 1/u<sub>t</sub>)*<sup>`kmers`</sup>), where *u<sub>t</sub>* is the number of specific *k*-mers for the tax id in the database.")
    public double getExpectedUniqueKMers() {
        return  (1 - Math.pow(1 - 1d / dbKMers, kmers)) * dbKMers;
    }

    @MDCDescription(pos = 16, name = "unique kmers / exp.", desc ="The ratio `unique kmers` / `exp. unique kmers` for the tax id. This should be close to 1 for a consistent match of *k*-mers. ([This paper](https://arxiv.org/pdf/1602.05822.pdf) discusses the corresponding background distribution (of `unique kmers`).)")
    public double getKMerConsistency() {
        return uniqueKmers / getExpectedUniqueKMers();
    }

    @MDCDescription(pos = 20, name ="db kmers", desc = "The number *u<sub>t</sub>* of specific *k*-mers for the tax id in the database.")
    public long getDbKMers() {
        return dbKMers;
    }

    @MDCDescription(pos = 21, name = "parent taxid", desc = "The parent tax id.")
    public String getParentTaxId() {
        return parentTaxId;
    }

    @MDCDescription(pos = 22, name = "mean error", desc = "The avg. ratio of inconsistent *k*-mers per read length for reads classified with respect to the tax id.")
    public double getMeanError() {
        return errorSum / reads;
    }

    @MDCDescription(pos = 23, name = "error std. dev.", desc = "The standard deviation of the `mean error`.")
    public double getErrorStdDev() {
        return Math.sqrt((errorSquaredSum  - errorSum * errorSum / reads) / (reads - 1));
    }

    public long getValueFor(ValueType valueType) {
        switch (valueType) {
            case READS:
                return reads;
            case READS_BPS:
                return readsBPs;
            case READS_1KMER:
                return reads1KMer;
            case READS_KMERS:
                return readsKmers;
            case KMERS:
                return kmers;
            default:
                throw new IllegalArgumentException();
        }
    }

    @MDCDescription(pos = 998, name="norm.", desc = "Normalized value for a respective value type.")
    public double getNormalizedFor(ValueType valueType) {
        return ((double) getValueFor(valueType)) / dbKMers;
    }

    @MDCDescription(pos = 999, name="acc.", desc = "Accumulated value or accumulated normalized valued for a respective value type.")
    public AccValues getAccValuesFor(ValueType valueType) {
        return extendedValues[valueType.ordinal()];
    }

    @MDCDescription(pos = 1000, name="max contig desc.", desc ="The descriptor of a read that holds a contiguous sequence of maximum length (according to the previous column).")
    public byte[] getMaxContigDescriptor() {
        return maxContigDescriptor;
    }

    @MDCDescription(pos = 1001, name = "acc. mean error", desc = "The accumulated `mean error`.")
    public double getAccMeanError() {
        CountsPerTaxid.AccValues accValues = getAccValuesFor(ValueType.READS);
        return accErrorSum / (accValues == null ? 0 : accValues.getAccumulated());
    }

    @MDCDescription(pos = 1002, name = "acc. error  std. dev.", desc = "The standard deviation of the `acc. mean error`.")
    public double getAccErrorDev() {
        CountsPerTaxid.AccValues accValues = getAccValuesFor(ValueType.READS);
        long reads = accValues == null ? 0 : accValues.getAccumulated();
        return Math.sqrt((accErrorSquaredSum  - (accErrorSum * accErrorSum) / reads) / (reads - 1));
    }

    @MDCDescription(pos = 2001, name="max kmer counts", desc = "The frequencies of the most frequent unique *k*-mers which are specific to the tax id's genome in descending order separated by `;`. " +
            "This column is experimental and only added when the configuration property `matchWithKMerCounts` is set to `true`. " +
            "The number of frequencies is determined via `maxKMerResCounts` (see also Section [Configuration parameters](#configuration-parameters)).")
    public short[] getMaxKMerCounts() {
        return maxKMerCounts;
    }

    void completeValues(int pos, long dbKMers, SmallTaxTree.SmallTaxIdNode node) {
        this.pos = pos;
        this.dbKMers = dbKMers;
        if (node != null) {
            this.name = node.getName();
            this.rank = node.getRank();
            this.parentTaxId = node.getParent() != null ? node.getParent().getTaxId() : "";
            for (int i = 0; i < ValueType.VALUES.length; i++) {
                long value = getValueFor(ValueType.VALUES[i]);
                extendedValues[i] = new AccValues(value, dbKMers);
            }
            accErrorSum = errorSum;
            accErrorSquaredSum = errorSquaredSum;
        }
        else {
            this.name = "TOTAL";
        }
    }

    void accumulateFrom(CountsPerTaxid other) {
        for (int i = 0; i < ValueType.VALUES.length; i++) {
            extendedValues[i].accumulateFrom(other.extendedValues[i]);
        }
        accErrorSum += other.accErrorSum;
        accErrorSquaredSum += other.accErrorSquaredSum;
    }
}