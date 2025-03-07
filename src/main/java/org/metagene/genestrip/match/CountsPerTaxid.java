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

import java.io.Serializable;

public class CountsPerTaxid implements Serializable {
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

    public static class  ExtendedValues implements Serializable {
        private final double normalized;
        private long accumulated;
        private double accumulatedNormalized;

        public ExtendedValues(long value, long dbKMers) {
            this.normalized = dbKMers > 0 ? ((double) value) / dbKMers : 0;
            this.accumulated = value;
            this.accumulatedNormalized = normalized;
        }

        public double getAccumulatedNormalized() {
            return accumulatedNormalized;
        }

        public double getNormalized() {
            return normalized;
        }

        public long getAccumulated() {
            return accumulated;
        }

        public void accumulateFrom(ExtendedValues extendedValues) {
            accumulated += extendedValues.getAccumulated();
            accumulatedNormalized += extendedValues.getAccumulatedNormalized();
        }
    };

    private static final long serialVersionUID = 1L;

    protected String taxid;
    protected long reads; // Only classified if a read is ENTIRELY consistent with the kmers of a taxon (and/or its ancestors) in the DB
    protected long reads1Kmer; // Reads with at least one k-mer belonging to this taxid
    protected long readsBPs;
    protected long readKmers; // Kmers from classified reads.
    protected long uniqueKmers; // All unique kmers counted - even from unclassified reads.
    protected long kmers; // All kmers counted - even from unclassified reads.
    protected int maxContigLen;
    protected int contigs;
    protected short[] maxKMerCounts;
    protected byte[] maxContigDescriptor;

    // For extended values
    private long dbKMers;
    private String parentTaxId;

    private final ExtendedValues[] extendedValues;

    public CountsPerTaxid(String taxid, int maxReadSizeBytes) {
        this.taxid = taxid;
        maxContigDescriptor = new byte[maxReadSizeBytes];
        extendedValues = new ExtendedValues[ValueType.VALUES.length];
    }

    public void setParentTaxId(String parentTaxId) {
        this.parentTaxId = parentTaxId;
    }

    public String getParentTaxId() {
        return parentTaxId;
    }

    public String getTaxid() {
        return taxid;
    }

    public int getContigs() {
        return contigs;
    }

    public long getKMers() {
        return kmers;
    }

    public int getMaxContigLen() {
        return maxContigLen;
    }

    public byte[] getMaxContigDescriptor() {
        return maxContigDescriptor;
    }

    public long getReads() {
        return reads;
    }

    public long getReadKMers() {
        return readKmers;
    }
    public long getUniqueKMers() {
        return uniqueKmers;
    }

    public short[] getMaxKMerCounts() {
        return maxKMerCounts;
    }

    public long getReads1Kmer() {
        return reads1Kmer;
    }

    public long getReadsBPs() {
        return readsBPs;
    }

    public double getAverageReadLength() {
        return (double) reads / readsBPs;
    }

    public double getCoverage() {
        return ((double) uniqueKmers) / dbKMers;
    }

    public double getExpectedUniqueKMers() {
        return  (1 - Math.pow(1 - 1d / dbKMers, kmers)) * dbKMers;
    }

    public long getDbKMers() {
        return dbKMers;
    }

    public long getForValueType(ValueType valueType) {
        switch (valueType) {
            case READS:
                return reads;
            case READS_BPS:
                return readsBPs;
            case READS_1KMER:
                return reads1Kmer;
            case READS_KMERS:
                return readKmers;
            case KMERS:
                return kmers;
            default:
                throw new IllegalArgumentException();
        }
    }

    void initExtendedValues(long dbKMers) {
        this.dbKMers = dbKMers;
        for (int i = 0; i < ValueType.VALUES.length; i++) {
            long value = getForValueType(ValueType.VALUES[i]);
            extendedValues[i] = new ExtendedValues(value, dbKMers);
        }
    }

    void accumulateFrom(CountsPerTaxid other) {
        for (int i = 0; i < ValueType.VALUES.length; i++) {
            extendedValues[i].accumulateFrom(other.extendedValues[i]);
        }
    }
}