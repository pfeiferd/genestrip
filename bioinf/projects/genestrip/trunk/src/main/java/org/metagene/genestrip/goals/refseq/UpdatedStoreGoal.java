package org.metagene.genestrip.goals.refseq;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;

public class UpdatedStoreGoal extends FilledStoreGoal {
	@SafeVarargs
	public UpdatedStoreGoal(GSProject project, String name, FileGoal<GSProject> updateStoreGoal,
			Goal<GSProject>... dependencies) {
		super(project, name, updateStoreGoal, Goal.append(dependencies, updateStoreGoal));
	}
}
