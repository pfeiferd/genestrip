package org.metagene.genestrip.make;

import java.io.File;
import java.io.IOException;
import java.util.List;

public abstract class FileGoal extends Goal {
	private final boolean allowEmptyFiles;

	public FileGoal(String name, Goal... dependencies) {
		this(name, false, dependencies);
	}

	public FileGoal(String name, boolean allowEmptyFiles, Goal... dependencies) {
		super(name, dependencies);
		this.allowEmptyFiles = allowEmptyFiles;
	}

	protected abstract List<File> getFiles();

	public boolean isAllowEmptyFiles() {
		return allowEmptyFiles;
	}

	@Override
	public boolean isMade() {
		for (File file : getFiles()) {
			if (!isMade(file)) {
				return false;
			}
		}
		return true;
	}

	protected boolean isMade(File file) {
		return file.exists() && (isAllowEmptyFiles() || file.length() != 0);
	}

	@Override
	public void makeThis() {
		try {
			for (File file : getFiles()) {
				if (!isMade(file)) {
					if (!isAllowEmptyFiles()) {
						if (file.exists() && file.length() == 0) {
							if (getLogger().isInfoEnabled()) {
								getLogger().info("Deleting emtpy file " + file);
							}
							file.delete();
						}
					}
					makeFile(file);
				}
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected abstract void makeFile(File file) throws IOException;

	@Override
	public void cleanThis() {
		for (File file : getFiles()) {
			if (file.exists()) {
				if (getLogger().isInfoEnabled()) {
					getLogger().info("Deleting file " + file);
				}
				file.delete();
			}
		}
	}

	@Override
	public String toString() {
		return "file goal: "+ getName();
	}
}
