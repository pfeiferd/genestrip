package org.metagene.genestrip.make;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public abstract class FileListGoal extends FileGoal {
	private final List<File> files;
	
	public FileListGoal(String name, Goal... dependencies) {
		this(name, null, false, dependencies);
	}
	
	public FileListGoal(String name, File file, Goal... dependencies) {
		this(name, Collections.singletonList(file), false, dependencies);
	}
	
	public FileListGoal(String name, List<File> files, Goal... dependencies) {
		this(name, files, false, dependencies);
	}
	
	public FileListGoal(String name, List<File> files, boolean allowEmptyFiles, Goal... dependencies) {
		super(name, dependencies);
		this.files = files != null ? new ArrayList<File>(files) : new ArrayList<File>();
	}
	
	public void addFile(File file) {
		files.add(file);
	}
	
	protected List<File> getFiles() {
		return files;
	}	
}
