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
    private static final long serialVersionUID = 1L;

    protected String taxid;
    protected long reads; // Only classified if a read is ENTIRELY consistent with the kmers of a taxon (and/or its ancestors) in the DB
    protected long reads1Kmer; // Reads with at least one k-mer belonging to this taxid
    protected long readsKmerBPs;
    protected long readKmers; // Kmers from classified reads.
    protected long uniqueKmers; // All unique kmers counted - even from unclassified reads.
    protected long kmers; // All kmers counted - even from unclassified reads.
    protected int maxContigLen;
    protected int contigs;
    protected short[] maxKMerCounts;
    protected byte[] maxContigDescriptor;

    private double noramlizedReads;
    private double accNoramlizedReads;
    private double coverage;
    private double expUnique;
    private long dbKMers;
    private long accReads;

    public CountsPerTaxid(String taxid, int maxReadSizeBytes) {
        this.taxid = taxid;
        maxContigDescriptor = new byte[maxReadSizeBytes];
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

    public long getReadsKmerBPs() {
        return readsKmerBPs;
    }

    public double getAccNormalizedReads() {
        return accNoramlizedReads;
    }

    public double getNormalizedReads() {
        return noramlizedReads;
    }

    public void addAccNormalizedReads(double accNoramlizedReads) {
        this.accNoramlizedReads += accNoramlizedReads;
    }

    public void setNormalizedReads(double noramlizedReads) {
        this.noramlizedReads = noramlizedReads;
    }

    public void setCoverage(double coverage) {
        this.coverage = coverage;
    }

    public double getCoverage() {
        return coverage;
    }

    public void setExpUnique(double expUnique) {
        this.expUnique = expUnique;
    }

    public double getExpUnique() {
        return expUnique;
    }

    public void setDbKMers(long dbKMers) {
        this.dbKMers = dbKMers;
    }

    public long getDbKMers() {
        return dbKMers;
    }

    public void addAccReads(long accReads) {
        this.accReads += accReads;
    }

    public long getAccReads() {
        return accReads;
    }
}