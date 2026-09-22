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
package org.metagene.genestrip.refseq;

import java.io.IOException;
import java.util.Collection;
import java.util.List;

import me.tongfei.progressbar.ProgressBar;
import org.apache.commons.logging.Log;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSConfigKey.SeqType;
import org.metagene.genestrip.io.BufferedLineReader;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.io.StreamingResource.StreamAccess;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.GSLogFactory;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;

/**
 * Streams and parses the (large) NCBI RefSeq accession catalog file line by line, filtering entries
 * by sequence type (DNA/RNA/mRNA), RefSeq category and release status, and dispatching each matching
 * entry to {@link #handleEntry(byte[], int, int, int)}.
 */
public abstract class AccessionFileProcessor {
    private static final int MAX_LINE_SIZE = 2048;

    /** Accession prefixes that mark any genomic (DNA) sequence, as against RNA and protein ones. */
    protected static final String[] GENOMIC_ACCESSION_PREFIXES = {"AC_", "NC_", "NG_", "NT_", "NW_", "NZ_"};

    // The subset of the above that a genome assembly is made of, which is a different question from
    // finished against draft: NZ_ marks every non-curated genomic sequence of the release, whole-
    // genome shotgun contigs included. What it leaves out are NG_, NT_ and NW_, the region records --
    // a few kilobases out of a genome, never an assembly of one. Hence two uses: whether a record's
    // sequence is wanted at all (`refseq.assemblyAccessionsOnly'), and whether it may take a place
    // under `maxGenomesPerTaxid', which counts genomes and so counts these alone.
    /** Accession prefixes of the kinds a genome assembly consists of. */
    protected static final String[] ASSEMBLY_ACCESSION_PREFIXES = {"AC_", "NC_", "NZ_"};

    /** Accession prefixes that mark a (non-messenger) RNA sequence. */
    protected static final String[] RNA_PREFIXES = {"NR_", "XR_"};

    /** Accession prefixes that mark a messenger-RNA sequence. */
    protected static final String[] M_RNA_PREFIXES = {"NM_", "XM_"};

    /** Logger used for accession-catalog reading. */
    protected final Log logger = GSLogFactory.getLog("accreader");

    private final RefSeqCategory[] categories;
    private final GSConfigKey.RefSeqStatus[] statuses;
    private final boolean dna;
    private final boolean rna;
    private final boolean mrna;

    /**
     * Creates a processor that keeps catalog entries matching the given RefSeq categories, sequence
     * type and release statuses.
     *
     * @param categories the RefSeq categories whose entries are kept
     * @param seqType the sequence type (DNA/RNA/mRNA) to keep
     * @param statuses the RefSeq release statuses to keep
     */
    public AccessionFileProcessor(Collection<RefSeqCategory> categories, SeqType seqType, List<GSConfigKey.RefSeqStatus> statuses) {
        // Converting to array makes iterator below way more efficient (less/no object
        // allocations) -
        // found via optimizer ...
        this.categories = categories.toArray(new RefSeqCategory[categories.size()]);
        dna = SeqType.GENOMIC.equals(seqType) || SeqType.ALL.equals(seqType);
        rna = SeqType.RNA.equals(seqType) || SeqType.ALL.equals(seqType) || SeqType.ALL_RNA.equals(seqType);
        mrna = SeqType.M_RNA.equals(seqType) || SeqType.ALL.equals(seqType) || SeqType.ALL_RNA.equals(seqType);
        this.statuses = statuses.toArray(new GSConfigKey.RefSeqStatus[statuses.size()]);
    }

    /**
     * Streams and parses the given RefSeq catalog file, invoking {@link #handleEntry} for every
     * entry that matches the configured sequence type, category and status filters.
     *
     * @param catalogFile the RefSeq catalog file to stream and parse
     */
    public void processCatalog(StreamingResource catalogFile) {
        try (StreamAccess byteCountAccess = catalogFile.openStream()) {
            long totalCatSize = byteCountAccess.getSize();
            // The file is huge, apache csv reader would be too slow and burn too many
            // strings. Therefore, manual coding for parsing and processing.
            byte[] target = new byte[MAX_LINE_SIZE];
            int size;
            try (ProgressBar pb = isProgressBar() ?
                    GSProgressBarCreator.newGSProgressBar(getProgressBarTaskName(), byteCountAccess, null) : null) {
                try (BufferedLineReader reader = new BufferedLineReader(byteCountAccess.getInputStream())) {
                    while ((size = reader.nextLine(target)) > 0) {
                        if (size > target.length) {
                            // nextLine returns target.length + 1 when a line does not fit the buffer;
                            // using it as a scan bound below would read past the array.
                            throw new IllegalStateException("buffer is too small for a line in the accession catalog file");
                        }
                        int pos1 = ByteArrayUtil.indexOf(target, 0, size, '\t');
                        int pos2 = ByteArrayUtil.indexOf(target, pos1 + 1, size, '\t');
                        int pos3 = ByteArrayUtil.indexOf(target, pos2 + 1, size, '\t');
                        int pos4 = ByteArrayUtil.indexOf(target, pos3 + 1, size, '\t');
                        int pos5 = ByteArrayUtil.indexOf(target, pos4 + 1, size, '\t');
                        // The line terminator is written into the buffer and counted in size, so the
                        // last column would otherwise be handed on with the newline still attached
                        // and would not parse as a number. Every other column ends at a tab and is
                        // unaffected; this one ends at the end of the line.
                        int lineEnd = size;
                        while (lineEnd > 0 && (target[lineEnd - 1] == '\n' || target[lineEnd - 1] == '\r')) {
                            lineEnd--;
                        }
                        if ((dna && isGenomicAccession(target, pos2 + 1)) || (rna && isRNAAccession(target, pos2 + 1))
                                || (mrna && isMRNAAccession(target, pos2 + 1))) {
                            if (containsCategory(target, pos3 + 1, pos4, categories)) {
                                if (containsStatus(target, pos4 + 1, pos5, statuses)) {
                                    handleEntry(target, pos1, pos2 + 1, pos3, pos5 + 1, lineEnd);
                                }
                            }
                        }
                   }
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Returns whether a progress bar should be shown while processing the catalog.
     *
     * @return {@code true} if a progress bar should be shown
     */
    protected boolean isProgressBar() {
        return true;
    }

    /**
     * Returns the task name shown on the progress bar.
     *
     * @return the progress bar task name
     */
    protected String getProgressBarTaskName() {
        return ((GSLogFactory.GSLog) logger).getName();
    }

    /**
     * Handles a catalog entry that passed all filters, given the line buffer and the byte offsets
     * delimiting its tax id (ending at {@code taxIdEnd}) and its accession ({@code accessionStart}
     * inclusive to {@code accessionEnd} exclusive).
     *
     * @param target the line buffer holding the catalog entry
     * @param taxIdEnd the exclusive end offset of the tax id
     * @param accessionStart the inclusive start offset of the accession
     * @param accessionEnd the exclusive end offset of the accession
     */
    protected abstract void handleEntry(byte[] target, int taxIdEnd, int accessionStart, int accessionEnd);

    /**
     * Handles one entry, with the length column too.
     * <p>
     * The catalog's last column is the length of the sequence, which the parse loop has already walked
     * past to find the accession, so offering it costs nothing. It is given as offsets rather than a
     * parsed number because most callers do not want it: the default implementation drops it and calls
     * {@link #handleEntry(byte[], int, int, int)}, so a subclass that does not care is unaffected.
     *
     * @param target the buffer holding the line
     * @param taxIdEnd the end offset of the tax id
     * @param accessionStart the start offset of the accession
     * @param accessionEnd the end offset of the accession
     * @param lengthStart the start offset of the length column
     * @param lengthEnd the end offset of the length column, i.e. the end of the line, with any line
     *            terminator already excluded
     */
    protected void handleEntry(byte[] target, int taxIdEnd, int accessionStart, int accessionEnd,
                               int lengthStart, int lengthEnd) {
        handleEntry(target, taxIdEnd, accessionStart, accessionEnd);
    }

    /**
     * Parses the non-negative number in {@code seq[start, end)}, or -1 where it is not one.
     *
     * @param seq the buffer holding the number
     * @param start the start offset
     * @param end the end offset
     * @return the number, or -1
     */
    public static long parseLength(byte[] seq, int start, int end) {
        if (start >= end) {
            return -1;
        }
        long v = 0;
        for (int i = start; i < end; i++) {
            byte c = seq[i];
            if (c < '0' || c > '9') {
                return -1;
            }
            v = v * 10 + (c - '0');
        }
        return v;
    }

    /**
     * Returns whether the given byte range contains the directory name of any of the given
     * categories.
     *
     * @param outerArray the buffer to search
     * @param start the inclusive start offset of the range
     * @param end the exclusive end offset of the range
     * @param categories the categories whose directory names are searched for
     * @return whether any category's directory name occurs in the range
     */
    protected boolean containsCategory(byte[] outerArray, int start, int end, RefSeqCategory[] categories) {
        for (int i = 0; i < categories.length; i++) {
            if (ByteArrayUtil.indexOf(outerArray, start, end, categories[i].getDirectory()) != -1) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns whether the given byte range contains the name of any of the given release statuses.
     *
     * @param outerArray the buffer to search
     * @param start the inclusive start offset of the range
     * @param end the exclusive end offset of the range
     * @param status the release statuses whose names are searched for
     * @return whether any status name occurs in the range
     */
    protected boolean containsStatus(byte[] outerArray, int start, int end, GSConfigKey.RefSeqStatus[] status) {
        for (int i = 0; i < status.length; i++) {
            if (ByteArrayUtil.indexOf(outerArray, start, end, status[i].getName()) != -1) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns whether the accession starting at the given offset has a genomic (DNA) prefix.
     *
     * @param outerArray the buffer holding the accession
     * @param start the start offset of the accession
     * @return whether the accession has a genomic prefix
     */
    protected boolean isGenomicAccession(byte[] outerArray, int start) {
        for (int i = 0; i < GENOMIC_ACCESSION_PREFIXES.length; i++) {
            if (ByteArrayUtil.startsWith(outerArray, start, GENOMIC_ACCESSION_PREFIXES[i])) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns whether the accession starting at the given offset has a (non-messenger) RNA prefix.
     *
     * @param outerArray the buffer holding the accession
     * @param start the start offset of the accession
     * @return whether the accession has an RNA prefix
     */
    public static boolean isRNAAccession(byte[] outerArray, int start) {
        for (int i = 0; i < RNA_PREFIXES.length; i++) {
            if (ByteArrayUtil.startsWith(outerArray, start, RNA_PREFIXES[i])) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns whether the accession starting at the given offset has a messenger-RNA prefix.
     *
     * @param outerArray the buffer holding the accession
     * @param start the start offset of the accession
     * @return whether the accession has a messenger-RNA prefix
     */
    public static boolean isMRNAAccession(byte[] outerArray, int start) {
        for (int i = 0; i < M_RNA_PREFIXES.length; i++) {
            if (ByteArrayUtil.startsWith(outerArray, start, M_RNA_PREFIXES[i])) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns whether the accession starting at the given offset has a complete-genome prefix.
     *
     * @param outerArray the buffer holding the accession
     * @param start the start offset of the accession
     * @return whether the accession has a complete-genome prefix
     */
    public static boolean isAssemblyAccession(byte[] outerArray, int start) {
        String[] prefixes = AccessionFileProcessor.ASSEMBLY_ACCESSION_PREFIXES;
        for (int i = 0; i < prefixes.length; i++) {
            if (ByteArrayUtil.startsWith(outerArray, start, prefixes[i])) {
                return true;
            }
        }
        return false;
    }
}
