package org.metagene.genestrip.goals;

import java.io.File;
import java.io.IOException;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.bloom.AbstractKMerBloomIndex;
import org.metagene.genestrip.bloom.FastqBloomFilter;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.util.ArraysUtil;

public class FilterGoal extends FileListGoal<GSProject> {
	private final File fastq;
	private final boolean writeDump;
	private final BloomFilterFileGoal bloomFilterFileGoal;

	@SafeVarargs
	public FilterGoal(GSProject project, String name, File fastq, boolean writeDump,
			BloomFilterFileGoal bloomFilterFileGoal, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile("filtered", fastq, FileType.FASTQ),
				ArraysUtil.append(deps, bloomFilterFileGoal));
		this.fastq = fastq;
		this.bloomFilterFileGoal = bloomFilterFileGoal;
		this.writeDump = writeDump;
	}

	@Override
	protected void makeFile(File file) {
		try {
			File dumpFile = writeDump ? getProject().getOutputFile("dumped", fastq, FileType.FASTQ) : null;
			new FastqBloomFilter(AbstractKMerBloomIndex.load(bloomFilterFileGoal.getOutputFile()),
					getProject().getConfig().getMinPosCountFilter(), getProject().getConfig().getPosRatioFilter(),
					getProject().getConfig().getMaxReadSizeBytes()).runFilter(fastq, file, dumpFile);
		} catch (IOException | ClassNotFoundException e) {
			throw new RuntimeException(e);
		}
	}
}
