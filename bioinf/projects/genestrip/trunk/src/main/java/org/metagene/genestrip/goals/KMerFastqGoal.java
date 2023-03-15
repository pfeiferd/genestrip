package org.metagene.genestrip.goals;

import java.io.File;
import java.io.IOException;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.fastqgen.KMerFastqGenerator;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;

public class KMerFastqGoal extends FileListGoal<GSProject> {
	private long addedKmers;
	private final FastasSizeGoal fastasSizeGoal;
	private final FastaFileDownloadGoal fastaDownloadGoal;

	@SafeVarargs
	public KMerFastqGoal(GSProject project, FastasSizeGoal fastasSizeGoal, FastaFileDownloadGoal fastaDownloadGoal,
			Goal<GSProject>... dependencies) {
		super(project, "kmerfastqgen", project.getKmerFastqFile(), dependencies);
		this.fastasSizeGoal = fastasSizeGoal;
		this.fastaDownloadGoal = fastaDownloadGoal;
	}

	@Override
	protected void makeFile(File fastq) {
		try {
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File generate " + fastq);
			}
			KMerFastqGenerator generator = new KMerFastqGenerator(getProject().getkMserSize());
			addedKmers = generator.run(fastaDownloadGoal.getFiles(), fastq, getProject().getName(),
					fastasSizeGoal.get());
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Entered KMers: " + addedKmers);
				getLogger().info("File generated " + fastq);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public long getAddedKmers() {
		return addedKmers;
	}
}
