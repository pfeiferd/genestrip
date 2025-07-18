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

import me.tongfei.progressbar.DelegatingProgressBarConsumer;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarBuilder;
import me.tongfei.progressbar.ProgressBarStyle;
import org.apache.commons.io.FileUtils;
import org.metagene.genestrip.GSConfigKey;

public abstract class FileGoal<P extends Project> extends Goal<P> {
	private final boolean allowEmptyFiles;

	@SafeVarargs
	public FileGoal(P project, GoalKey key, Goal<P>... dependencies) {
		this(project, key, false, dependencies);
	}

	@SafeVarargs
	public FileGoal(P project, GoalKey key, boolean allowEmptyFiles, Goal<P>... dependencies) {
		super(project, key, dependencies);
		this.allowEmptyFiles = allowEmptyFiles;
	}

	public abstract List<File> getFiles();

	public File getFile() {
		List<File> files = getFiles();
		if (files.size() == 1) {
			return files.get(0);
		}

		throw new IllegalStateException("There is zero or more than one output file here: " + files);
	}

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

	protected boolean isMade(File file) {
		return file.exists();
	}

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

	protected ProgressBar createProgressBar(int max) {
		return null;
	}

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

	protected void deleteFile(File file) {
		file.delete();
	}
	
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
