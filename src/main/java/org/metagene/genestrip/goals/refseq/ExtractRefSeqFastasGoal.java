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
            File file = getProject().getOutputFile(getKey().getName(), GSProject.GSFileType.CSV, false);
            try (PrintStream ps = new PrintStream(StreamProvider.getOutputStreamForFile(file))) {
                //ps.println("refseq descr; taxid;");
                for (String key : descr2TaxId.keySet()) {
                    ps.print(key);
                    ps.print("\t");
//                    ps.print(";");
                    ps.print(descr2TaxId.get(key));
//                    ps.println(";");
                    ps.println();
                }
            }
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
                regionsPerTaxid);
    }

    protected class MyFastaReader extends AbstractRefSeqFastaReader {
        private OutputStream os;

        public MyFastaReader(int bufferSize, Set<TaxTree.TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
                             int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId, int stepSize, boolean completeGenomesOnly, StringLong2DigitTrie regionsPerTaxid) {
            super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId, stepSize, completeGenomesOnly, regionsPerTaxid);
        }

        @Override
        protected void infoLine() {
            super.infoLine();
            if (includeRegion) {
                int pos = ByteArrayUtil.indexOf(target, 0, size, ' ');
                String name = new String(target, 1, pos - 1);
                String taxid = node.getTaxId();
                descr2TaxId.put(name, taxid);
                File file = new File(getProject().getFastaDir(), name + ".fasta.gz");
                try {
                    os = StreamProvider.getOutputStreamForFile(file);
                    PrintStream ps = new PrintStream(os);
                    ps.print('>');
                    ps.print(name);
                    ps.print(" ");
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
