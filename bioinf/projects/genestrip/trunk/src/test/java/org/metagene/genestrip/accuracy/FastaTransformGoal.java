package org.metagene.genestrip.accuracy;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.util.ByteArrayUtil;

public class FastaTransformGoal extends FileListGoal<GSProject> {

	private static final String[] NAMES = { "A_hydrophila_HiSeq", "B_cereus_HiSeq", "B_fragilis_HiSeq",
			"M_abscessus_HiSeq", "P_fermentans_HiSeq", "R_sphaeroides_HiSeq", "S_aureus_HiSeq", "S_pneumoniae_HiSeq",
			"V_cholerae_HiSeq", "X_axonopodis_HiSeq" };
	private static final String[] TAXIDS = { "1073377", "1053231", "1073387", "1001740", "1149860", "272943", "1213734",
			"170187", "991923", "1185664" };

	@SafeVarargs
	public FastaTransformGoal(GSProject project, String name, Goal<GSProject>... dependencies) {
		super(project, name, (List<File>) null, dependencies);
		addFile(new File(project.getFastqDir(), "HiSeq_accuracy.fastq"));
//		addFile(new File(project.getFastqDir(), "MiSeq_accuracy.fastq"));
	}

	@Override
	protected void makeFile(File file) throws IOException {
		PrintStream out = new PrintStream(new FileOutputStream(file));

		AbstractFastaReader fastaReader = new AbstractFastaReader(getProject().getConfig().getMaxReadSizeBytes()) {
			private int counter = 0;
			private int dataSize;

			@Override
			protected void infoLine() throws IOException {
				if (counter > 0) {
					printAdditionalLines();
				}
				dataSize = 0;
				int index = getNameIndexFromInfoLine();
				if (index != -1) {
					out.print("@");
					out.print(TAXIDS[index]);
					out.print(":");
					out.println(counter++);
				} else {
					if (getLogger().isWarnEnabled()) {
						getLogger().warn("Insconsistent info line: " + new String(target, 0, size));
					}
				}
			}

			@Override
			protected void dataLine() throws IOException {
				ByteArrayUtil.print(target, 0, size - 1, out);
				dataSize += size - 1;
			}

			@Override
			protected void done() throws IOException {
				printAdditionalLines();
			}

			protected void printAdditionalLines() {
				out.println();
				out.println("+");
				for (int i = 0; i < dataSize; i++) {
					out.print('~');
				}
				out.println();
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
		int dotIndex = file.getName().indexOf('.');
		String prefix = file.getName().substring(0, dotIndex);
		fastaReader.readFasta(new File(getProject().getFastaDir(), "accuracy/" + prefix + ".fa"));
	}
}
