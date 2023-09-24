package org.metagene.genestrip.goals.refseq;

import java.io.IOException;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerStoreWrapper;

public class FilledStoreGoal extends ObjectGoal<KMerStoreWrapper, GSProject> {
	private final FileGoal<GSProject> fillStoreGoal;

	@SafeVarargs
	public FilledStoreGoal(GSProject project, String name, FileGoal<GSProject> fillStoreGoal,
			Goal<GSProject>... dependencies) {
		super(project, name, Goal.append(dependencies, fillStoreGoal));
		this.fillStoreGoal = fillStoreGoal;
	}

	@Override
	public void makeThis() {
		KMerStoreWrapper wrapper;
		try {
			wrapper = KMerStoreWrapper.load(fillStoreGoal.getFile());
			set(wrapper);
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}		
	}
	
	void setStoreWrapper(KMerStoreWrapper object) {
		set(object);
	}
}
