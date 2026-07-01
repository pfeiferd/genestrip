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

/**
 * Holds all match statistics accumulated for a single tax id during read matching, such
 * as the number of classified reads and base pairs, matched and unique k-mer counts,
 * contig statistics and error measures. Getter methods annotated with
 * {@link MDCDescription} define the columns of the CSV match result. Instances are
 * ordered by their sort position within the taxonomy tree.
 */
public class CountsPerTaxid implements Serializable, Comparable<CountsPerTaxid> {
    /**
     * The kinds of countable values that are tracked and can be normalized and accumulated.
     */
    public enum ValueType {
        /** The number of classified reads. */
        READS("reads"),
        /** The number of matched k-mers. */
        KMERS("kmers"),
        /** The total number of base pairs of classified reads. */
        READS_BPS("reads bps"),
        /** The number of reads with at least one matched k-mer. */
        READS_1KMER("read >=1 kmer"),
        /** The number of k-mers from classified reads. */
        READS_KMERS("reads kmers");

        private static final ValueType[] VALUES = ValueType.values();

        private final String name;

        ValueType(String name) {
            this.name = name;
        }

        /**
         * Returns the display name of this value type.
         *
         * @return the display name
         */
        public String getName() {
            return name;
        }
    }

    /**
     * A raw value together with its version normalized by the number of database k-mers,
     * accumulated over a subtree of the taxonomy.
     */
    public static class AccValues implements Serializable {
        /** The accumulated raw value. */
        private long accumulated;
        /** The accumulated value normalized by the number of database k-mers. */
        private double accumulatedNormalized;

        /**
         * Creates accumulated values from a single raw value, normalizing by
         * {@code dbKMers} (or 0 if {@code dbKMers} is not positive).
         *
         * @param value   the raw value
         * @param dbKMers the number of database k-mers used for normalization
         */
        public AccValues(long value, long dbKMers) {
            this.accumulated = value;
            this.accumulatedNormalized = dbKMers > 0 ? ((double) value) / dbKMers : 0;
        }

        /**
         * Returns the accumulated normalized value.
         *
         * @return the accumulated normalized value
         */
        public double getAccumulatedNormalized() {
            return accumulatedNormalized;
        }

        /**
         * Returns the accumulated raw value.
         *
         * @return the accumulated raw value
         */
        public long getAccumulated() {
            return accumulated;
        }

        /**
         * Adds the raw and normalized values of the given instance into this one.
         *
         * @param extendedValues the values to add into this instance
         */
        public void accumulateFrom(AccValues extendedValues) {
            accumulated += extendedValues.getAccumulated();
            accumulatedNormalized += extendedValues.getAccumulatedNormalized();
        }
    }

    private static final long serialVersionUID = 1L;

    /** The depth of the tax id in the taxonomy tree. */
    protected int level;
    /** The tax id these counts belong to. */
    protected String taxid;
    /** The number of reads classified to the tax id. */
    protected long reads;
    /** The number of reads with at least one k-mer of the tax id. */
    protected long reads1KMer;
    /** The total number of base pairs of classified reads. */
    protected long readsBPs;
    /** The number of k-mers from classified reads consistent with the tax id. */
    protected long readsKmers;
    /** The number of unique matched k-mers (including those from unclassified reads). */
    protected long uniqueKmers; // All unique kmers counted - even from unclassified reads.
    /** The number of matched k-mers specific to the tax id's genome. */
    protected long kmers;
    /** The number of contiguous k-mer sequences specific to the tax id's genome. */
    protected int contigs;
    /** The sum of squared contig lengths, used for the contig-length standard deviation. */
    protected long contigLenSquaredSum;
    /** The maximum contig length in k-mers. */
    protected int maxContigLen;
    /** The read descriptor of the longest contig. */
    protected final byte[] maxContigDescriptor;
    /** The frequencies of the most frequent unique k-mers in descending order. */
    protected short[] maxKMerCounts;
    /** The sum of per-read k-mer error ratios. */
    protected double errorSum;
    /** The sum of squared per-read k-mer error ratios. */
    protected double errorSquaredSum;
    /** The sum of per-read class error ratios. */
    protected double classErrorSum;
    /** The sum of squared per-read class error ratios. */
    protected double classErrorSquaredSum;

    // To complete values
    /** The sort position of this entry within the taxonomy tree. */
    private int pos;
    /** The name associated with the tax id. */
    private String name;
    /** The taxonomic rank of the tax id. */
    private Rank rank;
    /** The number of specific k-mers for the tax id in the database. */
    private long dbKMers;
    /** The parent tax id. */
    private String parentTaxId;
    /** The accumulated values per value type. */
    private final AccValues[] extendedValues;
    /** The accumulated sum of per-read k-mer error ratios. */
    protected double accErrorSum;
    /** The accumulated sum of squared per-read k-mer error ratios. */
    protected double accErrorSquaredSum;
    /** The accumulated sum of per-read class error ratios. */
    protected double accClassErrorSum;
    /** The accumulated sum of squared per-read class error ratios. */
    protected double accClassErrorSquaredSum;

    /**
     * Creates an empty entry for the given tax id.
     *
     * @param level            the depth of the tax id in the taxonomy tree
     * @param taxid            the tax id
     * @param maxReadSizeBytes the maximum read descriptor size in bytes
     */
    public CountsPerTaxid(int level, String taxid, int maxReadSizeBytes) {
        this.level = level;
        this.taxid = taxid;
        maxContigDescriptor = new byte[maxReadSizeBytes];
        extendedValues = new AccValues[ValueType.VALUES.length];
    }

    /**
     * Creates an entry pre-populated with total counts.
     *
     * @param level         the depth of the tax id in the taxonomy tree
     * @param taxid         the tax id
     * @param totalReads    the total number of reads
     * @param totalKmers    the total number of k-mers
     * @param totalBPs      the total number of base pairs
     * @param totalMaxCounts the maximum k-mer counts
     */
    public CountsPerTaxid(int level, String taxid, long totalReads, long totalKmers, long totalBPs, short[] totalMaxCounts) {
        this.level = level;
        this.taxid = taxid;
        this.reads = totalReads;
        this.kmers = totalKmers;
        this.readsBPs = totalBPs;
        this.maxKMerCounts = totalMaxCounts;
        maxContigDescriptor = new byte[0];
        extendedValues = new AccValues[ValueType.VALUES.length];
    }

    /**
     * Orders entries by their sort position within the taxonomy tree.
     */
    @Override
    public int compareTo(CountsPerTaxid o) {
        return this.pos - o.pos;
    }

    /**
     * Returns the depth in the taxonomy tree, where the root is zero.
     *
     * @return the depth in the taxonomy tree
     */
    @MDCDescription(pos = 0, name = "level", desc = "Depth in the taxonomy tree, where the root is zero.")
    public int getLevel() {
        return level;
    }

    /**
     * Returns the sort position of this entry.
     *
     * @return the sort position
     */
    @MDCDescription(pos = -1, name = "pos", desc = "Sort position of entry.")
    public int getPos() {
        return pos;
    }

    /**
     * Returns the name associated with the tax id.
     *
     * @return the name associated with the tax id
     */
    @MDCDescription(pos = 1, name = "name", desc = "The name associated with the tax id.")
    public String getName() {
        return name;
    }

    /**
     * Returns the rank of the tax id.
     *
     * @return the rank of the tax id
     */
    @MDCDescription(pos = 2, name= "rank", desc = "The rank of the tax id.")
    public Rank getRank() {
        return rank;
    }

    /**
     * Returns the tax id.
     *
     * @return the tax id
     */
    @MDCDescription(pos = 3, name = "taxid", desc = "The tax id.")
    public String getTaxid() {
        return taxid;
    }

    // Only classified if a read is ENTIRELY consistent with the kmers of a taxon (and/or its ancestors) in the DB
    /**
     * Returns the number of reads classified to the tax id.
     *
     * @return the number of classified reads
     */
    @MDCDescription(pos = 4, name = "reads", desc = "The number of reads classified with respect to the tax id.")
    public long getReads() {
        return reads;
    }

    /**
     * Returns the number of k-mers from classified reads consistent with the read's tax id.
     *
     * @return the number of k-mers from classified reads
     */
    @MDCDescription(pos = 5, name="kmers from reads", desc ="The number of *k*-mers from classified reads which are consistent with the read's tax id.")
    public long getReadsKMers() {
        return readsKmers;
    }

    /**
     * Returns all matched k-mers specific to the tax id's genome.
     *
     * @return the number of matched k-mers
     */
    @MDCDescription(pos = 6, name="kmers", desc = "*All* matched *k*-mers which are specific to the tax id's genome (according to the database). The *k*-mers do not have to be in an accordingly classified read for this.")
    public long getKMers() {
        return kmers;
    }

    /**
     * Returns all unique matched k-mers specific to the tax id's genome.
     *
     * @return the number of unique matched k-mers
     */
    @MDCDescription(pos = 7, name="unique kmers", desc = "*All* unique *k*-mers, which are specific to the tax id's genome (according to the database). " +
            "Here, multiple occurrences of the same *k*-mer are only counted once. " +
            "Genestrip always performs exact counting according to [KrakenUniq's exact counting](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1568-0#Sec8). (Genestrip implements an efficient in-memory storage method for related counts based on a bit vector.)")
    public long getUniqueKMers() {
        return uniqueKmers;
    }

    /**
     * Returns the number of contiguous k-mer sequences specific to the tax id's genome.
     *
     * @return the number of contigs
     */
    @MDCDescription(pos = 8, name="contigs", desc = "The number of contiguous sequences of *k*-mers that are specific to the tax id's genome.")
    public int getContigs() {
        return contigs;
    }

    /**
     * Returns the average contig length in k-mers.
     *
     * @return the average contig length
     */
    @MDCDescription(pos = 9, name="average contig length", desc = "The average length of contiguous sequences of *k*-mers (as the number of *k*-mers) that are specific to the tax id's genome.")
    public double getAverageContigLen() {
        return ((double) kmers) / contigs;
    }

    /**
     * Returns the maximum contig length in k-mers.
     *
     * @return the maximum contig length
     */
    @MDCDescription(pos = 10, name="max contig length", desc ="The maximum length of all contiguous sequences of *k*-mers (as the number of *k*-mers) that are specific to the tax id's genome.")
    public int getMaxContigLen() {
        return maxContigLen;
    }

    /**
     * Returns the number of reads with at least one k-mer of the tax id.
     *
     * @return the number of reads with at least one k-mer
     */
    @MDCDescription(pos = 11, name = "reads >=1 kmer", desc = "Reads with at least one *k*-mer of the respective tax id.")
    public long getReads1KMer() {
        return reads1KMer;
    }

    /**
     * Returns the total number of base pairs of classified reads.
     *
     * @return the total number of base pairs
     */
    @MDCDescription(pos = 12, name = "reads bps", desc = "The total number of base pairs of reads classified with respect to the tax id.")
    public long getReadsBPs() {
        return readsBPs;
    }

    /**
     * Returns the average length of classified reads in base pairs.
     *
     * @return the average read length in base pairs
     */
    @MDCDescription(pos = 13, name = "avg. read length", desc = "The average length of classified reads in base pairs.")
    public double getAverageReadLength() {
        return ((double) readsBPs) / reads;
    }

    /**
     * Returns the database coverage as unique k-mers divided by db k-mers.
     *
     * @return the database coverage
     */
    @MDCDescription(pos = 14, name = "db coverage", desc ="The ratio `unique kmers` / u<sub>t</sub>, where *u<sub>t</sub>* = `db kmers`")
    public double getCoverage() {
        return ((double) uniqueKmers) / dbKMers;
    }

    /**
     * Returns the expected number of unique k-mers.
     *
     * @return the expected number of unique k-mers
     */
    @MDCDescription(pos = 15, name = "exp. unique kmers", desc ="The number of expected unique *k*-mers, which is *u<sub>t</sub> * (1 - (1 - 1/u<sub>t</sub>)*<sup>`kmers`</sup>), where *u<sub>t</sub>* is the number of specific *k*-mers for the tax id in the database.")
    public double getExpectedUniqueKMers() {
        return  (1 - Math.pow(1 - 1d / dbKMers, kmers)) * dbKMers;
    }

    /**
     * Returns the ratio of unique k-mers to expected unique k-mers.
     *
     * @return the k-mer consistency ratio
     */
    @MDCDescription(pos = 16, name = "unique kmers / exp.", desc ="The ratio `unique kmers` / `exp. unique kmers` for the tax id. This should be close to 1 for a consistent match of *k*-mers. ([This paper](https://arxiv.org/pdf/1602.05822.pdf) discusses the corresponding background distribution (of `unique kmers`).)")
    public double getKMerConsistency() {
        return uniqueKmers / getExpectedUniqueKMers();
    }

    /**
     * Returns the number of specific k-mers for the tax id in the database.
     *
     * @return the number of database k-mers
     */
    @MDCDescription(pos = 20, name ="db kmers", desc = "The number *u<sub>t</sub>* of specific *k*-mers for the tax id in the database.")
    public long getDbKMers() {
        return dbKMers;
    }

    /**
     * Returns the parent tax id.
     *
     * @return the parent tax id
     */
    @MDCDescription(pos = 21, name = "parent taxid", desc = "The parent tax id.")
    public String getParentTaxId() {
        return parentTaxId;
    }

    /**
     * Returns the mean k-mer error over classified reads.
     *
     * @return the mean k-mer error
     */
    @MDCDescription(pos = 22, name = "mean error", desc = "The mean ratio of a classified read's *k*-mers that are not in the database per read's total *k*-mers.")
    public double getMeanError() {
        return errorSum / reads;
    }

    /**
     * Returns the standard deviation of the mean error.
     *
     * @return the standard deviation of the mean error
     */
    @MDCDescription(pos = 23, name = "kmer error std. dev.", desc = "The standard deviation of the `mean error`.")
    public double getErrorStdDev() {
        return Math.sqrt((errorSquaredSum  - errorSum * errorSum / reads) / (reads - 1));
    }

    /**
     * Returns the mean class error over reads.
     *
     * @return the mean class error
     */
    @MDCDescription(pos = 24, name = "mean class error", desc = "The mean ratio of a read's *k*-mers that are not consistent with the read's class per read's total *k*-mers.")
    public double getMeanClassError() {
        return classErrorSum / reads;
    }

    /**
     * Returns the standard deviation of the mean class error.
     *
     * @return the standard deviation of the mean class error
     */
    @MDCDescription(pos = 25, name = "class error std. dev.", desc = "The standard deviation of the `mean class error`.")
    public double getClassErrorStdDev() {
        return Math.sqrt((classErrorSquaredSum  - classErrorSum * classErrorSum / reads) / (reads - 1));
    }

    /**
     * Returns the standard deviation of the contig length.
     *
     * @return the standard deviation of the contig length
     */
    @MDCDescription(pos = 26, name = "contig len std. dev.", desc = "The standard deviation of the `contig length`.")
    public double getContigLenStdDev() {
        return Math.sqrt((contigLenSquaredSum  - ((double) kmers * kmers) / contigs) / (contigs - 1));
    }

    /**
     * Returns the raw count for the given value type.
     *
     * @param valueType the value type whose raw count is requested
     * @return the raw count for the value type
     * @throws IllegalArgumentException if the value type is not a raw countable value
     */
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

    /**
     * Returns the normalized value for the given value type.
     *
     * @param valueType the value type whose normalized value is requested
     * @return the normalized value
     */
    @MDCDescription(pos = 998, name="norm.", desc = "Normalized value of a respective value type.")
    public double getNormalizedFor(ValueType valueType) {
        return ((double) getValueFor(valueType)) / dbKMers;
    }

    /**
     * Returns the accumulated values for the given value type.
     *
     * @param valueType the value type whose accumulated values are requested
     * @return the accumulated values
     */
    @MDCDescription(pos = 999, name="acc.", desc = "Accumulated value or accumulated normalized valued of a respective value type.")
    public AccValues getAccValuesFor(ValueType valueType) {
        return extendedValues[valueType.ordinal()];
    }

    /**
     * Returns the descriptor of the read holding the longest contig.
     *
     * @return the descriptor of the longest contig
     */
    @MDCDescription(pos = 1000, name="max contig desc.", desc ="The descriptor of a read that holds a contiguous sequence of maximum length (according to the column `max contig length`).")
    public byte[] getMaxContigDescriptor() {
        return maxContigDescriptor;
    }

    /**
     * Returns the accumulated mean error.
     *
     * @return the accumulated mean error
     */
    @MDCDescription(pos = 1001, name = "acc. mean error", desc = "The accumulated `mean error`.")
    public double getAccMeanError() {
        CountsPerTaxid.AccValues accValues = getAccValuesFor(ValueType.READS);
        return accErrorSum / (accValues == null ? 0 : accValues.getAccumulated());
    }

    /**
     * Returns the standard deviation of the accumulated mean error.
     *
     * @return the standard deviation of the accumulated mean error
     */
    @MDCDescription(pos = 1002, name = "acc. error std. dev.", desc = "The standard deviation of the `acc. mean error`.")
    public double getAccErrorStdDev() {
        CountsPerTaxid.AccValues accValues = getAccValuesFor(ValueType.READS);
        long reads = accValues == null ? 0 : accValues.getAccumulated();
        return Math.sqrt((accErrorSquaredSum  - (accErrorSum * accErrorSum) / reads) / (reads - 1));
    }

    /**
     * Returns the accumulated mean class error.
     *
     * @return the accumulated mean class error
     */
    @MDCDescription(pos = 1003, name = "acc. mean class error", desc = "The accumulated `mean class error`.")
    public double getAccClassMeanError() {
        CountsPerTaxid.AccValues accValues = getAccValuesFor(ValueType.READS);
        return accClassErrorSum / (accValues == null ? 0 : accValues.getAccumulated());
    }

    /**
     * Returns the standard deviation of the accumulated mean class error.
     *
     * @return the standard deviation of the accumulated mean class error
     */
    @MDCDescription(pos = 1004, name = "acc. class error std. dev.", desc = "The standard deviation of the `acc. mean class error`.")
    public double getAccClassStdErrorDev() {
        CountsPerTaxid.AccValues accValues = getAccValuesFor(ValueType.READS);
        long reads = accValues == null ? 0 : accValues.getAccumulated();
        return Math.sqrt((accClassErrorSquaredSum  - (accClassErrorSum * accClassErrorSum) / reads) / (reads - 1));
    }

    /**
     * Returns the frequencies of the most frequent unique k-mers in descending order.
     *
     * @return the maximum k-mer counts
     */
    @MDCDescription(pos = 2001, name="max kmer counts", desc = "The frequencies of the most frequent unique *k*-mers which are specific to the tax id's genome in descending order separated by `;`. " +
            "This column is experimental and only present when the configuration parameter `maxKMerResCounts` is set to a value greater than 0 " +
            "(see also Section [Configuration parameters](#configuration-parameters)).")
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
            accClassErrorSum = classErrorSum;
            accClassErrorSquaredSum = classErrorSquaredSum;
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
        accClassErrorSum += other.accClassErrorSum;
        accClassErrorSquaredSum += other.accClassErrorSquaredSum;
    }
}