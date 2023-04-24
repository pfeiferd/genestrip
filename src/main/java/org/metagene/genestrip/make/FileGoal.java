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
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.FileUtils;

public abstract class FileGoal<P> extends Goal<P> {
	private final boolean allowEmptyFiles;

	@SafeVarargs
	public FileGoal(P project, String name, Goal<P>... dependencies) {
		this(project, name, false, dependencies);
	}

	@SafeVarargs
	public FileGoal(P project, String name, boolean allowEmptyFiles, Goal<P>... dependencies) {
		super(project, name, dependencies);
		this.allowEmptyFiles = allowEmptyFiles;
	}

	protected abstract List<File> getFiles();

	public File getOutputFile() {
		List<File> files = getFiles();
		if (getFiles().size() == 1) {
			return files.get(0);
		}
		return null;
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

	protected boolean isMade(File file) {
		return file.exists();
	}
	
	protected boolean isBadEmpty(File file) {
		return !isAllowEmptyFiles() && !file.isDirectory() && file.length() == 0;
	}

	@Override
	public void makeThis() {
		try {
			for (File file : getFiles()) {
				if (isBadEmpty(file)) {
					if (getLogger().isInfoEnabled()) {
						getLogger().info("Deleting emtpy file " + file);
					}
					file.delete();
				}
				if (!isMade(file)) {
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
		try {
			for (File file : getFilesToClean()) {
				if (file.exists()) {
					if (getLogger().isInfoEnabled()) {
						getLogger().info("Deleting file " + file);
					}
					if (file.isDirectory()) {
						FileUtils.deleteDirectory(file);
					} else {
						file.delete();
					}
				}
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	@Override
	protected boolean isThisCleaned() {
		for (File file : getFilesToClean()) {
			if (file.exists()) {
				return false;
			}
		}
		return true;
	}

	protected List<File> getFilesToClean() {
		return getFiles();
	}

	@Override
	public String toString() {
		return "file goal: " + getName();
	}
}
