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
package org.metagene.genestrip.goals;

import me.tongfei.progressbar.ProgressBar;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.make.*;
import org.metagene.genestrip.match.ResultReporter;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.store.InvalidDatabaseClassException;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Properties;

public class DBConfigInfoGoal<P extends GSProject> extends ObjectGoal<Properties, P> {
	private final ObjectGoal<Database, P> dbGoal;
	private final File dbFile;

	@SafeVarargs
	public DBConfigInfoGoal(P project, GoalKey key, ObjectGoal<Database, P> dbGoal, FileGoal<P> storedDBGoal, Goal<P>... deps) {
		super(project, key,
				project.getDBPath() == null ? append(deps, dbGoal, storedDBGoal) : deps);
		this.dbGoal = dbGoal;
		this.dbFile = getProject().getDBPath() == null ? storedDBGoal.getFile() : new File(getProject().getDBPath());
	}

	@Override
	protected void doMakeThis() {
		try {
			if (dbGoal.isMade()) {
				set(dbGoal.get().getConfigInfo());
			} else {
				set(Database.loadConfigInfo(dbFile));
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
