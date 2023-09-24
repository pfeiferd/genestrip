package org.metagene.genestrip.goals.refseq;

import java.io.IOException;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;

public class BloomIndexedGoal extends ObjectGoal<MurmurCGATBloomFilter, GSProject> {
	private final FileGoal<GSProject> bloomIndexGoal;

	@SafeVarargs
	public BloomIndexedGoal(GSProject project, String name, FileGoal<GSProject> bloomIndexGoal,
			Goal<GSProject>... dependencies) {
		super(project, name, Goal.append(dependencies, bloomIndexGoal));
		this.bloomIndexGoal = bloomIndexGoal;
	}

	@Override
	public void makeThis() {
		try {
			MurmurCGATBloomFilter filter = MurmurCGATBloomFilter.load(bloomIndexGoal.getFile());
			set(filter);
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}		
	}
	
	void setFilter(MurmurCGATBloomFilter object) {
		set(object);
	}
}
