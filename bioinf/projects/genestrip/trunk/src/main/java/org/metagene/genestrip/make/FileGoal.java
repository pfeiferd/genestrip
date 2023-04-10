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
		for (File file : getFilesToClean()) {
			if (file.exists()) {
				if (getLogger().isInfoEnabled()) {
					getLogger().info("Deleting file " + file);
				}
				file.delete();
			}
		}
	}
	
	protected List<File> getFilesToClean() {
		return getFiles();
	}

	@Override
	public String toString() {
		return "file goal: "+ getName();
	}
}
