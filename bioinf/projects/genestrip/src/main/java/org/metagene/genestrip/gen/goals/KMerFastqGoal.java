package org.metagene.genestrip.gen.goals;

import java.io.File;
import java.io.IOException;

import org.metagene.genestrip.fastqgen.KMerFastqGenerator;
import org.metagene.genestrip.gen.Project;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;

public class KMerFastqGoal extends FileListGoal {
	private final Project project;
	private long addedKmers;
	private final FastasSizeGoal fastasSizeGoal;
	private final FastaFileDownloadGoal fastaDownloadGoal;

	public KMerFastqGoal(Project project, FastasSizeGoal fastasSizeGoal,
			FastaFileDownloadGoal fastaDownloadGoal, Goal... dependencies) {
		super("kmerfastqgen", project.getKmerFastqFile(), dependencies);
		this.project = project;
		this.fastasSizeGoal = fastasSizeGoal;
		this.fastaDownloadGoal = fastaDownloadGoal;
	}

	@Override
	protected void makeFile(File fastq) {
		try {
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File generate " + fastq);
			}
			KMerFastqGenerator generator = new KMerFastqGenerator(project.getkMserSize());
			addedKmers = generator.run(fastaDownloadGoal.getFiles(), fastq, project.getName(), fastasSizeGoal.get());
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
