package org.metagene.genestrip.gen.goals;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.metagene.genestrip.fastqgen.KrakenKMerFastqMerger;
import org.metagene.genestrip.fastqgen.KrakenKMerFastqMerger.FilterListener;
import org.metagene.genestrip.gen.Project;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class KrakenFastqFileGoal extends FileListGoal {
	private final Project project;
	private final ObjectGoal<Set<TaxIdNode>> taxNodesGoal;

	public KrakenFastqFileGoal(Project project, ObjectGoal<Set<TaxIdNode>> taxNodesGoal, Goal... deps) {
		super("krakenfastq", project.getFilteredKmerFastqFile(), deps);
		this.project = project;
		this.taxNodesGoal = taxNodesGoal;
	}

	@Override
	protected void makeFile(File fastq) {
		try {
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Writing file " + fastq);
			}
			PrintStream filteredOut = new PrintStream(new GZIPOutputStream(new FileOutputStream(fastq), 4096));

			FilterListener filter = KrakenKMerFastqMerger.createFilterByTaxIdNodes(taxNodesGoal.get(),
					KrakenKMerFastqMerger.createFastQOutputFilterByTaxId(filteredOut, null));
			KrakenKMerFastqMerger krakenKMerFastqMerger = new KrakenKMerFastqMerger(
					project.getConfig().getMaxReadSizeBytes());

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Reading file " + project.getKrakenOutFile());
				getLogger().info("Reading file " + project.getKmerFastqFile());
			}
			FileInputStream fStream = new FileInputStream(project.getKmerFastqFile());
			GZIPInputStream gStream = new GZIPInputStream(fStream, 4096);
			krakenKMerFastqMerger.process(new BufferedInputStream(new FileInputStream(project.getKrakenOutFile())),
					new BufferedInputStream(gStream), filter);

			filteredOut.close();
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Written file " + project.getFilteredKmerFastqFile());
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

}
