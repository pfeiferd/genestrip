package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.util.Collections;
import java.util.List;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;

public class RefSeqRNumDownloadGoal extends RefSeqDownloadGoal {
	public static final String RELEASE_NUMBER_FILE_NAME = "RELEASE_NUMBER";
	
	private final List<File> files;
	
	@SafeVarargs
	public RefSeqRNumDownloadGoal(GSProject project, String name, Goal<GSProject>... deps) {
		super(project, name, deps);
		
		files = Collections.singletonList(new File(project.getConfig().getRefSeqDir(), RELEASE_NUMBER_FILE_NAME));
	}

	@Override
	protected String getFTPDir(File file) {
		return RELEASE_FOLDER;
	}
	
	@Override
	public List<File> getFiles() {
		return files;
	}
}
