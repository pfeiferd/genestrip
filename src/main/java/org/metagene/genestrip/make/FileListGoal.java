package org.metagene.genestrip.make;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public abstract class FileListGoal<P> extends FileGoal<P> {
	private final List<File> files;
	
	@SafeVarargs
	public FileListGoal(P project, String name, Goal<P>... dependencies) {
		this(project, name, null, false, dependencies);
	}
	
	@SafeVarargs
	public FileListGoal(P project, String name, File file, Goal<P>... dependencies) {
		this(project, name, Collections.singletonList(file), false, dependencies);
	}
	
	@SafeVarargs
	public FileListGoal(P project, String name, List<File> files, Goal<P>... dependencies) {
		this(project, name, files, false, dependencies);
	}
	
	@SafeVarargs
	public FileListGoal(P project, String name, List<File> files, boolean allowEmptyFiles, Goal<P>... dependencies) {
		super(project, name, dependencies);
		this.files = files != null ? new ArrayList<File>(files) : new ArrayList<File>();
	}
	
	public void addFile(File file) {
		files.add(file);
	}
	
	protected List<File> getFiles() {
		return files;
	}	
}
