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

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.fastqgen.KrakenKMerFastqMerger;
import org.metagene.genestrip.fastqgen.KrakenKMerFastqMerger.FilterListener;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class KrakenFastqFileGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;

	@SafeVarargs
	public KrakenFastqFileGoal(GSProject project, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			Goal<GSProject>... deps) {
		super(project, "krakenfastq", project.getFilteredKmerFastqFile(), deps);
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
					getProject().getConfig().getMaxReadSizeBytes());

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Reading file " + getProject().getKrakenOutFile());
				getLogger().info("Reading file " + getProject().getKmerFastqFile());
			}
			FileInputStream fStream = new FileInputStream(getProject().getKmerFastqFile());
			GZIPInputStream gStream = new GZIPInputStream(fStream, 4096);
			krakenKMerFastqMerger.process(new BufferedInputStream(new FileInputStream(getProject().getKrakenOutFile())),
					new BufferedInputStream(gStream), filter);

			filteredOut.close();
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Written file " + getProject().getFilteredKmerFastqFile());
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

}
