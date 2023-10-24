package org.metagene.genestrip.accuracy;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.util.ByteArrayUtil;

public class FastaTransformGoal extends FileListGoal<GSProject> {
	private static final String[] NAMES = { "A_hydrophila_HiSeq", "B_cereus_HiSeq", "B_fragilis_HiSeq",
			"M_abscessus_HiSeq", "P_fermentans_HiSeq", "R_sphaeroides_HiSeq", "S_aureus_HiSeq", "S_pneumoniae_HiSeq",
			"V_cholerae_HiSeq", "X_axonopodis_HiSeq",
			//
			"B_cereus_MiSeq", "M_abscessus_MiSeq", "R_sphaeroides_MiSeq", "V_cholerae_MiSeq", "C_freundii_MiSeq",
			"E_cloacae_MiSeq", "K_pneumoniae_MiSeq", "P_vulgaris_MiSeq", "S_aureus_MiSeq", "S_enterica_MiSeq" };
	public static final String[] TAXIDS = { "1073377", "1053231", "1073387", "1001740", "1149860", "272943", "1213734",
			"170187", "991923", "1185664",
			//
			"1053231", "1001740", "272943", "991923", "1173724", "550", "1173763", "1173773", "1210042", "1173880" };

	private final ObjectGoal<Map<String, String>, GSProject> mapGoal;

	@SafeVarargs
	public FastaTransformGoal(GSProject project, String name, ObjectGoal<Map<String, String>, GSProject> mapGoal,
			Goal<GSProject>... dependencies) {
		super(project, name, (List<File>) null, append(dependencies, mapGoal));
		this.mapGoal = mapGoal;
		addFile(new File(project.getFastqDir(), "HiSeq_accuracy.fastq"));
		addFile(new File(project.getFastqDir(), "MiSeq_accuracy.fastq"));
		addFile(new File(project.getFastqDir(), "simBA5_accuracy.fastq"));
		addFile(new File(project.getFastqDir(), "HiSeq_timing.fastq.gz"));
		addFile(new File(project.getFastqDir(), "MiSeq_timing.fastq.gz"));
		addFile(new File(project.getFastqDir(), "simBA5_timing.fastq.gz"));
	}

	@Override
	protected void makeFile(File file) throws IOException {
		int dotIndex = file.getName().indexOf('.');
		String prefix = file.getName().substring(0, dotIndex);
		File fastaFile = new File(getProject().getFastaDir(),
				(prefix.contains("accuracy") ? "accuracy/" : "") + prefix + ".fa");

		PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));

		boolean accNumberWay = prefix.startsWith("simBA5");

		final Map<String, String> accesion2TaxidMap = accNumberWay ? mapGoal.get() : null;

		AbstractFastaReader fastaReader = new AbstractFastaReader(getProject().getConfig().getMaxReadSizeBytes()) {
			private int counter = 0;
			private int dataSize;
			private boolean taxidFound;

			@Override
			protected void infoLine() throws IOException {
				if (taxidFound) {
					printAdditionalLines();
				}
				taxidFound = false;
				dataSize = 0;
				if (accNumberWay) {
					String accNumber = AccessionNumber2TaxidGoal.getAccessionNumberFromInfoLine(target, size);
					if (accNumber != null) {
						String taxId = accesion2TaxidMap.get(accNumber);
						if (taxId != null) {
							printDescriptor(taxId);
						} else {
							if (getLogger().isWarnEnabled()) {
								getLogger().warn("Missing taxid for accession number: " + accNumber);
							}
						}
					} else {
						if (getLogger().isWarnEnabled()) {
							getLogger().warn("Inconsistent info line: " + new String(target, 0, size));
						}
					}
				} else {
					int index = getNameIndexFromInfoLine();
					if (index != -1) {
						printDescriptor(TAXIDS[index]);
					} else {
						if (getLogger().isWarnEnabled()) {
							getLogger().warn("Inconsistent info line: " + new String(target, 0, size));
						}
					}
				}
			}

			@Override
			protected void dataLine() throws IOException {
				if (taxidFound) {
					ByteArrayUtil.print(target, 0, size - 1, out);
				}
				dataSize += size - 1;
			}

			@Override
			protected void done() throws IOException {
				printAdditionalLines();
			}

			protected void printDescriptor(String taxid) {
				taxidFound = true;
				out.print("@");
				out.print(taxid);
				out.print(":");
				out.println(counter++);
			}

			protected void printAdditionalLines() {
				if (taxidFound) {
					out.println();
					out.println("+");
					for (int i = 0; i < dataSize; i++) {
						out.print('~');
					}
					out.println();
				}
			}

			protected int getNameIndexFromInfoLine() {
				for (int i = 0; i < NAMES.length; i++) {
					if (ByteArrayUtil.startsWith(target, 1, NAMES[i])) {
						return i;
					}
				}
				return -1;
			}

		};

		fastaReader.readFasta(fastaFile);
	}
}
