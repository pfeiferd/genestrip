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
package org.metagene.genestrip.finertree.goals;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.ResultReporter;
import org.metagene.genestrip.store.Database;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

/**
 * Writes a CSV file containing general information about the database, as produced by
 * {@link ResultReporter#printStoreInfo}.
 *
 * @param <P> the concrete FT project type
 */
public class FTDBInfoGoal<P extends FTProject> extends FileListGoal<P> {
	private final ObjectGoal<Database, P> storeGoal;

	/**
	 * Creates the goal writing the database-info CSV file.
	 *
	 * @param project   the FT project
	 * @param storeGoal the goal providing the database
	 * @param deps      further goals this goal depends on
	 */
	@SafeVarargs
	public FTDBInfoGoal(P project, ObjectGoal<Database, P> storeGoal, Goal<P>... deps) {
		super(project, FTGoalKey.FTDBINFO, getDBInfoFile(project),
				Goal.append(deps, storeGoal));
		this.storeGoal = storeGoal;
	}

	/**
	 * Returns the database-info CSV output file for the given project.
	 *
	 * @param project the project
	 * @param <P>     the concrete project type
	 * @return the database-info CSV output file
	 */
	public static <P extends GSProject> File getDBInfoFile(P project) {
		return project.getOutputFile(FTGoalKey.FTDBINFO.getName(), GSProject.GSFileType.CSV, false);
	}

	/**
	 * Writes general information about the database to the given file.
	 *
	 * @param file the file to write
	 */
	@Override
	protected void makeFile(File file) {
		try {
			Database wrapper = storeGoal.get();
			try (PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file))) {
				new ResultReporter().printStoreInfo(wrapper, out,
						booleanConfigValue(GSConfigKey.DB_INFO_WITH_DISTANCE));
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
