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

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.*;

public class ExtractRefSeqFastasGoal<P extends GSProject> extends FastaReaderGoal<Map<String, String>, P> {
    private final ObjectGoal<AccessionMap, P> accessionMapGoal;
    private final List<FillSizeGoal.MyFastaReader> readers;

    private Map<String, String> descr2TaxId;

    @SafeVarargs
    public ExtractRefSeqFastasGoal(P project, ExecutionContext bundle, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
                                   ObjectGoal<Set<TaxTree.TaxIdNode>, P> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
                                   ObjectGoal<AccessionMap, P> accessionMapGoal, Goal<P>... deps) {
        super(project, GSGoalKey.EXTRACT_REFSEQ_FASTA, bundle, categoriesGoal, taxNodesGoal, fnaFilesGoal, null, Goal.append(deps, accessionMapGoal));
        this.accessionMapGoal = accessionMapGoal;
        readers = new ArrayList<>();
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
            cleanUpThreads();
        }
    }

    @Override
    protected AbstractRefSeqFastaReader createFastaReader(AbstractRefSeqFastaReader.StringLong2DigitTrie regionsPerTaxid) {
        return new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
                taxNodesGoal.get(), isIncludeRefSeqFna() ? accessionMapGoal.get() : null, intConfigValue(GSConfigKey.KMER_SIZE),
                intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
                (Rank) configValue(GSConfigKey.MAX_GENOMES_PER_TAXID_RANK),
                longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID),
                intConfigValue(GSConfigKey.STEP_SIZE),
                booleanConfigValue(GSConfigKey.COMPLETE_GENOMES_ONLY),
                regionsPerTaxid,
                booleanConfigValue(GSConfigKey.EXTRACT_REFSEQ_GZIP));
    }

    protected class MyFastaReader extends AbstractRefSeqFastaReader {
        private OutputStream os;
        private final boolean gzip;

        public MyFastaReader(int bufferSize, Set<TaxTree.TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
                             int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid, boolean gzip) {
            super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, stepSize, completeGenomesOnly, regionsPerTaxid);
            this.gzip = gzip;
        }

        @Override
        protected void infoLine() {
            super.infoLine();
            if (includeRegion) {
                int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
                String name = new String(target, 1, pos - 1);
                String taxid = node.getTaxId();
                descr2TaxId.put(name, taxid);
                File file = new File(getProject().getFastaDir(), name + (gzip ? ".fa.gz" : ".fa"));
                try {
                    os = StreamProvider.getOutputStreamForFile(file);
                    PrintStream ps = new PrintStream(os);
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
            if (includeRegion) {
                try {
                    os.write(target, 0, size);
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }

        protected void endRegion() {
            if (includeRegion) {
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
