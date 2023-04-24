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
	private final BloomFilterFileGoal bloomFilterFileGoal;

	@SafeVarargs
	public FilterGoal(GSProject project, String name, BloomFilterFileGoal bloomFilterFileGoal,
			Goal<GSProject>... deps) {
		super(project, name, project.getOutputFileForFastq("filtered", FileType.FASTQ),
				ArraysUtil.append(deps, bloomFilterFileGoal));
		this.bloomFilterFileGoal = bloomFilterFileGoal;
	}

	@Override
	protected void makeFile(File file) {
		try {
			new FastqBloomFilter(AbstractKMerBloomIndex.load(bloomFilterFileGoal.getOutputFile()),
					getProject().getConfig().getMinPosCountFilter(), getProject().getConfig().getPosRatioFilter(),
					getProject().getConfig().getMaxReadSizeBytes()).runFilter(getProject().getFastqFile(), file,
							getProject().getDumpFastqFile());
		} catch (IOException | ClassNotFoundException e) {
			throw new RuntimeException(e);
		}
	}
}
