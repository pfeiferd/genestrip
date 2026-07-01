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

/**
 * A {@link FileGoal} whose set of output files is an explicit, mutable list. The files may be
 * supplied at construction time or added lazily via {@link #provideFiles()} the first time they are
 * requested.
 *
 * @param <P> the type of {@link Project} this goal belongs to
 */
public abstract class FileListGoal<P extends Project> extends FileGoal<P> {
	private boolean filesProvided;
	private final List<File> files;

	/**
	 * Creates a file list goal with a single output file.
	 *
	 * @param project      the project this goal belongs to
	 * @param key          the key identifying this goal
	 * @param file         the single output file produced by this goal
	 * @param dependencies the goals this goal depends on
	 */
	@SafeVarargs
	public FileListGoal(P project, GoalKey key, File file, Goal<P>... dependencies) {
		this(project, key, Collections.singletonList(file), false, dependencies);
	}

	/**
	 * Creates a file list goal with the given output files.
	 *
	 * @param project      the project this goal belongs to
	 * @param key          the key identifying this goal
	 * @param files        the output files produced by this goal
	 * @param dependencies the goals this goal depends on
	 */
	@SafeVarargs
	public FileListGoal(P project, GoalKey key, List<File> files, Goal<P>... dependencies) {
		this(project, key, files, false, dependencies);
	}

	/**
	 * Creates a file list goal with the given output files and empty-file policy. A {@code null}
	 * list starts an empty list to be populated later via {@link #addFile(File)}.
	 *
	 * @param project         the project this goal belongs to
	 * @param key             the key identifying this goal
	 * @param files           the initial output files, or {@code null} to start with an empty list
	 * @param allowEmptyFiles whether empty output files are permitted
	 * @param dependencies    the goals this goal depends on
	 */
	@SafeVarargs
	public FileListGoal(P project, GoalKey key, List<File> files, boolean allowEmptyFiles, Goal<P>... dependencies) {
		super(project, key, allowEmptyFiles, dependencies);
		this.files = files != null ? new ArrayList<File>(files) : new ArrayList<File>();
	}

	/**
	 * Adds a file to the output file list; typically called from {@link #provideFiles()}.
	 *
	 * @param file the output file to add
	 */
	protected void addFile(File file) {
		files.add(file);
	}

	/**
	 * Hook to lazily add output files the first time {@link #getFiles()} is called; empty by default.
	 */
	protected void provideFiles() {
	}

	/**
	 * Returns the output files, invoking {@link #provideFiles()} once on first access.
	 */
	public List<File> getFiles() {
		if (!filesProvided) {
			provideFiles();
			filesProvided = true;
		}
		return Collections.unmodifiableList(files);
	}
}
