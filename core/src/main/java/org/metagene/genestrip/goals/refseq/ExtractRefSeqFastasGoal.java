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
package org.metagene.genestrip.goals.refseq;
import java.nio.charset.StandardCharsets;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.tax.TaxNodeSelection;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.*;

/**
 * Goal that extracts the selected contigs into individual FASTA files (with kraken2-style
 * {@code |kraken:taxid|} headers) under the project's FASTA directory, and produces a map from each
 * sequence description to its taxid.
 * <p>
 * The selection matches the one of {@link FillDBGoal}: the same RefSeq categories, the same
 * additional FASTA files -- those downloaded from Genbank and those named in the project's
 * {@code additional.txt} -- and the same per-tax-id limits on contigs and $k$-mers. The extracted
 * files therefore represent exactly the genomes a database is filled from, which is what makes them
 * a valid basis for simulating reads against that database.
 *
 * @param <P> the project type
 */
public class ExtractRefSeqFastasGoal<P extends GSProject> extends FastaReaderGoal<Map<String, String>, P> {
    private final ObjectGoal<AccessionMap, P> accessionMapGoal;
    private Map<String, String> descr2TaxId;

    /**
     * Creates the goal, wiring the categories, tax-node, RefSeq-file, additional-FASTA and
     * accession-map goals it reads.
     *
     * @param project the project type
     * @param bundle the execution context providing threading and shared services
     * @param categoriesGoal the goal supplying the selected RefSeq categories
     * @param taxNodesGoal the goal supplying the selected taxonomic nodes
     * @param fnaFilesGoal the goal supplying the downloaded RefSeq FASTA files
     * @param additionalGoal the goal supplying additional FASTA files mapped to their tax node,
     *                       i.e. the Genbank downloads and the project's own additional files
     * @param accessionMapGoal the goal supplying the accession-to-tax-id map
     * @param deps the additional goals this goal depends on
     */
    @SafeVarargs
    public ExtractRefSeqFastasGoal(P project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                                   ObjectGoal<TaxNodeSelection, P> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                                   ObjectGoal<Map<File, TaxTree.TaxIdNode>, P> additionalGoal,
                                   ObjectGoal<AccessionMap, P> accessionMapGoal, Goal<P>... deps) {
        super(project, GSGoalKey.EXTRACT_REFSEQ_FASTA, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, additionalGoal, Goal.append(deps, accessionMapGoal));
        this.accessionMapGoal = accessionMapGoal;
        descr2TaxId = Collections.synchronizedMap(new HashMap<>());
    }

    @Override
    protected void doMakeThis() {
        try {
            readFastas();
            set(descr2TaxId);
        } catch (IOException e) {
            throw new RuntimeException(e);
        } finally {
            descr2TaxId = null;
        }
    }

    @Override
    protected AbstractRefSeqFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie contigsPerTaxid) {
        return new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
                taxNodesGoal.get().getSelected(), isIncludeRefSeqFna() ? accessionMapGoal.get() : null, intConfigValue(GSConfigKey.KMER_SIZE),
                intConfigValue(GSConfigKey.KMER_SAMPLING),
                booleanConfigValue(GSConfigKey.ASSEMBLY_ACCESSIONS_ONLY),
                contigsPerTaxid,
                booleanConfigValue(GSConfigKey.EXTRACT_REFSEQ_GZIP));
    }

    /**
     * FASTA reader that writes each included contig to its own FASTA file and records the mapping
     * from sequence description to taxid.
     */
    protected class MyFastaReader extends AbstractRefSeqFastaReader {
        private OutputStream os;
        private final boolean gzip;

        /**
         * Creates the reader that writes each included contig to its own FASTA file.
         *
         * @param bufferSize the FASTA line read-buffer size in bytes
         * @param taxNodes the taxonomic nodes to keep contigs for
         * @param accessionMap the accession-to-tax-id map, or {@code null} if not used
         * @param k the k-mer size
         * @param kMerSampling the k-mer sampling step size
         * @param assemblyAccessionsOnly whether only genomic accessions are considered, dropping `NG_`, `NT_` and `NW_`
         * @param contigsPerTaxid the per-tax-id contig counter
         * @param gzip whether the output FASTA files are GZIP-compressed
         */
        public MyFastaReader(int bufferSize, Set<TaxTree.TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
                             int kMerSampling, boolean assemblyAccessionsOnly, StringLong2DigitTrie contigsPerTaxid, boolean gzip) {
            super(bufferSize, taxNodes, accessionMap, k, kMerSampling, assemblyAccessionsOnly, contigsPerTaxid);
            this.gzip = gzip;
        }

        @Override
        protected void infoLine() {
            super.infoLine();
            if (includeContig) {
                // The name is the accession, i.e. everything up to the description. A header without
                // a description has none of the latter, so the contig reaches to the end of the line.
                int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
                if (pos < 0) {
                    pos = size;
                    while (pos > 1 && (target[pos - 1] == '\n' || target[pos - 1] == '\r')) {
                        pos--;
                    }
                }
                String name = new String(target, 1, pos - 1, StandardCharsets.UTF_8);
                String taxid = node.getTaxId();
                descr2TaxId.put(name, taxid);
                File file = new File(getProject().getFastaDir(), name + (gzip ? ".fa.gz" : ".fa"));
                try {
                    os = StreamProvider.getOutputStreamForFile(file);
                    PrintStream ps = new PrintStream(os, false, StandardCharsets.UTF_8);
                    ps.print('>');
                    ps.print(name);
                    ps.print("|kraken:taxid|"); // This is to please kraken2 during library building...
                    ps.println(taxid);
                    ps.flush();
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }

        @Override
        protected void dataLine() {
            if (includeContig) {
                try {
                    os.write(target, 0, size);
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }

        protected void endContig() {
            if (includeContig) {
                try {
                    os.close();
                    os = null;
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }
    }
}
