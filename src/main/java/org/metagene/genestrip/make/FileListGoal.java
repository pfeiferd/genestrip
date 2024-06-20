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
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public abstract class FileListGoal<P extends Project> extends FileGoal<P> {
	private boolean filesProvided;
	private final List<File> files;

	@SafeVarargs
	public FileListGoal(P project, GoalKey key, File file, Goal<P>... dependencies) {
		this(project, key, Collections.singletonList(file), false, dependencies);
	}

	@SafeVarargs
	public FileListGoal(P project, GoalKey key, List<File> files, Goal<P>... dependencies) {
		this(project, key, files, false, dependencies);
	}

	@SafeVarargs
	public FileListGoal(P project, GoalKey key, List<File> files, boolean allowEmptyFiles, Goal<P>... dependencies) {
		super(project, key, dependencies);
		this.files = files != null ? new ArrayList<File>(files) : new ArrayList<File>();
	}

	protected void addFile(File file) {
		files.add(file);
	}

	protected void provideFiles() {
	}

	public List<File> getFiles() {
		if (!filesProvided) {
			provideFiles();
			filesProvided = true;
		}
		return Collections.unmodifiableList(files);
	}
}
