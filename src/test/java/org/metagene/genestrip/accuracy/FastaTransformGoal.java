package org.metagene.genestrip.accuracy;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.io.BufferedLineReader;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.util.ByteArrayUtil;

public class FastaTransformGoal extends FileListGoal<GSProject> {
	public static final String ACCESSION_MAP_FILE = "nucl_gb.accession2taxid.gz";
	public static final String OLD_ACCESSION_MAP_FILE = "dead_nucl.accession2taxid.gz";

	private static final String[] NAMES = { "A_hydrophila_HiSeq", "B_cereus_HiSeq", "B_fragilis_HiSeq",
			"M_abscessus_HiSeq", "P_fermentans_HiSeq", "R_sphaeroides_HiSeq", "S_aureus_HiSeq", "S_pneumoniae_HiSeq",
			"V_cholerae_HiSeq", "X_axonopodis_HiSeq",
			//
			"B_cereus_MiSeq", "M_abscessus_MiSeq", "R_sphaeroides_MiSeq", "V_cholerae_MiSeq", "C_freundii_MiSeq",
			"E_cloacae_MiSeq", "K_pneumoniae_MiSeq", "P_vulgaris_MiSeq", "S_aureus_MiSeq", "S_enterica_MiSeq" };
	private static final String[] TAXIDS = { "1073377", "1053231", "1073387", "1001740", "1149860", "272943", "1213734",
			"170187", "991923", "1185664",
			//
			"1053231", "1001740", "272943", "991923", "1173724", "550", "1173763", "1173773", "1210042", "1173880" };

	private static final String MARKER = "|ref|";

	@SafeVarargs
	public FastaTransformGoal(GSProject project, String name, Goal<GSProject>... dependencies) {
		super(project, name, (List<File>) null, dependencies);
		addFile(new File(project.getFastqDir(), "HiSeq_accuracy.fastq"));
		addFile(new File(project.getFastqDir(), "MiSeq_accuracy.fastq"));
		addFile(new File(project.getFastqDir(), "simBA5_accuracy.fastq"));
	}

	@Override
	protected void makeFile(File file) throws IOException {
		int dotIndex = file.getName().indexOf('.');
		String prefix = file.getName().substring(0, dotIndex);
		File fastaFile = new File(getProject().getFastaDir(), "accuracy/" + prefix + ".fa");

		PrintStream out = new PrintStream(new FileOutputStream(file));

		boolean accNumberWay = prefix.startsWith("simBA5");

		final Map<String, String> accesion2TaxidMap = new HashMap<String, String>();
		if (accNumberWay) {
			Set<String> accNumbers = accNumberWay ? new HashSet<String>() : null;
			AbstractFastaReader fastaReader1 = new AbstractFastaReader(getProject().getConfig().getMaxReadSizeBytes()) {
				@Override
				protected void infoLine() throws IOException {
					String acc = getAccessionNumberFromInfoLine(target, size);
					if (acc != null) {
						accNumbers.add(acc);
					}
				}

				@Override
				protected void dataLine() throws IOException {
				}
			};
			fastaReader1.readFasta(fastaFile);
			fillAccessionNumbersToTaxIdsMap(ACCESSION_MAP_FILE, accesion2TaxidMap, accNumbers);
			fillAccessionNumbersToTaxIdsMap(OLD_ACCESSION_MAP_FILE, accesion2TaxidMap, accNumbers);
		}

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
					String accNumber = getAccessionNumberFromInfoLine(target, size);
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

	protected String getAccessionNumberFromInfoLine(byte[] target, int size) {
		int pos1 = ByteArrayUtil.indexOf(target, 1, size, MARKER);
		if (pos1 != -1) {
			pos1 += MARKER.length();
			int pos2 = ByteArrayUtil.indexOf(target, pos1, size, '|');
			if (pos2 != -1) {
				return new String(target, pos1, pos2 - pos1).trim();
			}
		}
		return null;
	}

	protected void fillAccessionNumbersToTaxIdsMap(String fileName, Map<String, String> accesion2TaxidMap, Set<String> accNumbers)
			throws IOException {

		InputStream in = StreamProvider
				.getInputStreamForFile(new File(getProject().getConfig().getCommonDir(), fileName));

		BufferedLineReader br = new BufferedLineReader(in);
		int size;
		byte[] target = new byte[2048];
		br.nextLine(target);
		while ((size = br.nextLine(target)) > 0) {
			int tab1 = ByteArrayUtil.indexOf(target, 0, size, '\t');
			int tab2 = ByteArrayUtil.indexOf(target, tab1 + 1, size, '\t');
			String acc = new String(target, tab1 + 1, tab2 - tab1 - 1);
			if (accNumbers.contains(acc)) {
				int tab3 = ByteArrayUtil.indexOf(target, tab2 + 1, size, '\t');
				String tax = new String(target, tab2 + 1, tab3 - tab2 - 1);
				accesion2TaxidMap.put(acc, tax);
			}
		}
	}
}
