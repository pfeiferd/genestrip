/*
 * 
 * “Commons Clause” License Condition v1.0
 * 
 * The Software is provided to you by the Licensor under the License, 
 * as defined below, subject to the following condition.
 * 
 * Without limiting other conditions in the License, the grant of rights under the License 
 * will not include, and the License does not grant to you, the right to Sell the Software.
 * 
 * For purposes of the foregoing, “Sell” means practicing any or all of the rights granted 
 * to you under the License to provide to third parties, for a fee or other consideration 
 * (including without limitation fees for hosting or consulting/ support services related to 
 * the Software), a product or service whose value derives, entirely or substantially, from the 
 * functionality of the Software. Any license notice or attribution required by the License 
 * must also include this Commons Clause License Condition notice.
 * 
 * Software: genestrip
 * 
 * License: Apache 2.0
 * 
 * Licensor: Daniel Pfeifer (daniel.pfeifer@progotec.de)
 * 
 */
package org.metagene.genestrip.make;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import me.tongfei.progressbar.ProgressBar;
import org.apache.commons.io.FileUtils;

/**
 * A {@link Goal} whose output is one or more files on disk. It is considered made when all of its
 * output files exist (and, unless empty files are allowed, are non-empty), and it is cleaned by
 * deleting those files. Subclasses implement {@link #makeFile(File)} to produce each file.
 *
 * @param <P> the type of {@link Project} this goal belongs to
 */
public abstract class FileGoal<P extends Project> extends Goal<P> {
	private final boolean allowEmptyFiles;

	/**
	 * Creates a file goal that does not allow empty output files.
	 *
	 * @param project the project this goal belongs to
	 * @param key the key identifying this goal
	 * @param dependencies the goals that must be made before this goal
	 */
	@SafeVarargs
	public FileGoal(P project, GoalKey key, Goal<P>... dependencies) {
		this(project, key, false, dependencies);
	}

	/**
	 * Creates a file goal. If {@code allowEmptyFiles} is {@code false}, an existing zero-length
	 * output file is treated as not made and gets deleted and remade.
	 *
	 * @param project the project this goal belongs to
	 * @param key the key identifying this goal
	 * @param allowEmptyFiles whether zero-length output files are considered made
	 * @param dependencies the goals that must be made before this goal
	 */
	@SafeVarargs
	public FileGoal(P project, GoalKey key, boolean allowEmptyFiles, Goal<P>... dependencies) {
		super(project, key, dependencies);
		this.allowEmptyFiles = allowEmptyFiles;
	}

	/**
	 * Returns the list of output files produced by this goal.
	 *
	 * @return the output files produced by this goal
	 */
	public abstract List<File> getFiles();

	/**
	 * Returns the single output file of this goal.
	 *
	 * @return the single output file of this goal
	 * @throws IllegalStateException if this goal has zero or more than one output file
	 */
	public File getFile() {
		List<File> files = getFiles();
		if (files.size() == 1) {
			return files.get(0);
		}

		throw new IllegalStateException("There is zero or more than one output file here: " + files);
	}

	/**
	 * Returns whether zero-length output files are considered made for this goal.
	 *
	 * @return {@code true} if empty output files are allowed
	 */
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
	
	@Override
	protected void alreadyMade() {
		if (getLogger().isDebugEnabled()) {
			getLogger().debug("Already made files " + getFiles());
		}
	}

	/**
	 * Returns whether the given file is made, i.e. it exists and is not a disallowed empty file.
	 *
	 * @param file the output file to check
	 * @return {@code true} if the file is considered made
	 */
	protected boolean isMade(File file) {
		return file.exists() && !isBadEmpty(file);
	}

	/**
	 * Returns whether the given file is an empty file that is not allowed, i.e. a non-directory that
	 * exists and has zero length while empty files are disallowed.
	 *
	 * @param file the output file to check
	 * @return {@code true} if the file is a disallowed empty file
	 */
	protected boolean isBadEmpty(File file) {
		return !isAllowEmptyFiles() && !file.isDirectory() && file.exists() && file.length() == 0;
	}

	@Override
	protected void doMakeThis() {
		try {
			int counter = 0;
			List<File> files = getFiles();
			int max = files.size();
			try (ProgressBar pb = createProgressBar(max)) {
				for (File file : files) {
					counter++;
					if (isBadEmpty(file)) {
						if (getLogger().isInfoEnabled()) {
							getLogger().info("Deleting emtpy file " + file);
						}
						deleteFile(file);
					}
					if (!isMade(file)) {
						if (getLogger().isDebugEnabled()) {
							getLogger().debug("Making file (" + counter + "/" + max + "):" + file);
						}
						makeFile(file);
					} else {
						if (getLogger().isDebugEnabled()) {
							getLogger().debug("Already made file (" + counter + "/" + max + "):" + file);
						}
					}
					if (pb != null) {
						pb.step();
					}
				}
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Creates a progress bar for making the output files, or {@code null} for none (the default).
	 *
	 * @param max the number of files to be made
	 * @return the progress bar to use, or {@code null} for none
	 */
	protected ProgressBar createProgressBar(int max) {
		return null;
	}

	/**
	 * Produces the given output file.
	 *
	 * @param file the output file to produce
	 * @throws IOException if the file cannot be written
	 */
	protected abstract void makeFile(File file) throws IOException;

	@Override
	protected void doCleanThis() {
		try {
			for (File file : getFilesToClean()) {
				if (file.exists()) {
					if (getLogger().isInfoEnabled()) {
						getLogger().info("Deleting file " + file);
					}
					if (file.isDirectory()) {
						try {
							FileUtils.deleteDirectory(file);
						} catch (FileNotFoundException e) {
							// Ignore on purpose. Happens on my Mac for unknown detail reasons.
						}
					} else {
						deleteFile(file);
					}
				}
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Deletes the given file.
	 *
	 * @param file the file to delete
	 */
	protected void deleteFile(File file) {
		file.delete();
	}
	
	/**
	 * Returns the files to delete when cleaning this goal; by default the output files.
	 *
	 * @return the files to delete when cleaning this goal
	 */
	protected List<File> getFilesToClean() {
		return getFiles();
	}

	@Override
	public boolean isCleaned() {
		for (File file : getFilesToClean()) {
			if (file.exists()) {
				return false;
			}
		}
		return true;
	}

	@Override
	public String toString() {
		return "file goal: " + getKey().getName();
	}
}
