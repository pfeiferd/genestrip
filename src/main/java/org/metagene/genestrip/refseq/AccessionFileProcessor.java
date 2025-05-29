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
import java.io.InputStream;
import java.util.Collection;

import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarBuilder;
import me.tongfei.progressbar.ProgressBarStyle;
import org.apache.commons.logging.Log;
import org.metagene.genestrip.GSConfigKey.SeqType;
import org.metagene.genestrip.io.BufferedLineReader;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.io.StreamingResource.StreamAccess;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.GSLogFactory;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;

public abstract class AccessionFileProcessor {
    private static int MAX_LINE_SIZE = 2048;

    protected static final String[] ALL_GENOMIC_ACCESSION_PREFIXES = {"AC_", "NC_", "NG_", "NT_", "NW_", "NZ_"};

    protected static final String[] COMPLETE_GENOMIC_ACCESSION_PREFIXES = {"AC_", "NC_", "NZ_"};

    protected static final String[] RNA_PREFIXES = {"NR_", "XR_"};

    protected static final String[] M_RNA_PREFIXES = {"NM_", "XM_"};

    protected final Log logger = GSLogFactory.getLog("accreader");

    private final RefSeqCategory[] categories;
    private final boolean dna;
    private final boolean rna;
    private final boolean mrna;

    public AccessionFileProcessor(Collection<RefSeqCategory> categories, SeqType seqType) {
        // Converting to array makes iterator below way more efficient (less/no object
        // allocations) -
        // found via optimizer ...
        this.categories = categories.toArray(new RefSeqCategory[categories.size()]);
        dna = SeqType.GENOMIC.equals(seqType) || SeqType.ALL.equals(seqType);
        rna = SeqType.RNA.equals(seqType) || SeqType.ALL.equals(seqType) || SeqType.ALL_RNA.equals(seqType);
        mrna = SeqType.M_RNA.equals(seqType) || SeqType.ALL.equals(seqType) || SeqType.ALL_RNA.equals(seqType);
    }

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
                        int pos1 = ByteArrayUtil.indexOf(target, 0, size, '\t');
                        int pos2 = ByteArrayUtil.indexOf(target, pos1 + 1, size, '\t');
                        int pos3 = ByteArrayUtil.indexOf(target, pos2 + 1, size, '\t');
                        int pos4 = ByteArrayUtil.indexOf(target, pos3 + 1, size, '\t');
                        if ((dna && isGenomicAccession(target, pos2 + 1)) || (rna && isRNAAccession(target, pos2 + 1))
                                || (mrna && isMRNAAccession(target, pos2 + 1))) {
                            if (containsCategory(target, pos3 + 1, pos4, categories)) {
                                handleEntry(target, pos1, pos2 + 1, pos3);
                            }
                        }
                   }
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    protected boolean isProgressBar() {
        return true;
    }

    protected String getProgressBarTaskName() {
        return ((GSLogFactory.GSLog) logger).getName();
    }

    protected abstract void handleEntry(byte[] target, int taxIdEnd, int accessionStart, int accessionEnd);

    protected boolean containsCategory(byte[] outerArray, int start, int end, RefSeqCategory[] categories) {
        for (int i = 0; i < categories.length; i++) {
            if (ByteArrayUtil.indexOf(outerArray, start, end, categories[i].getDirectory()) != -1) {
                return true;
            }
        }
        return false;
    }

    protected boolean isGenomicAccession(byte[] outerArray, int start) {
        for (int i = 0; i < ALL_GENOMIC_ACCESSION_PREFIXES.length; i++) {
            if (ByteArrayUtil.startsWith(outerArray, start, ALL_GENOMIC_ACCESSION_PREFIXES[i])) {
                return true;
            }
        }
        return false;
    }

    public static boolean isRNAAccession(byte[] outerArray, int start) {
        for (int i = 0; i < RNA_PREFIXES.length; i++) {
            if (ByteArrayUtil.startsWith(outerArray, start, RNA_PREFIXES[i])) {
                return true;
            }
        }
        return false;
    }

    public static boolean isMRNAAccession(byte[] outerArray, int start) {
        for (int i = 0; i < M_RNA_PREFIXES.length; i++) {
            if (ByteArrayUtil.startsWith(outerArray, start, M_RNA_PREFIXES[i])) {
                return true;
            }
        }
        return false;
    }

    public static boolean isCompleteGenomicAccession(byte[] outerArray, int start) {
        String[] prefixes = AccessionFileProcessor.COMPLETE_GENOMIC_ACCESSION_PREFIXES;
        for (int i = 0; i < prefixes.length; i++) {
            if (ByteArrayUtil.startsWith(outerArray, start, prefixes[i])) {
                return true;
            }
        }
        return false;
    }
}
