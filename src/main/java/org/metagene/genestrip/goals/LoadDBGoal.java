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

import java.io.File;
import java.io.IOException;
import java.io.InvalidClassException;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;

public class LoadDBGoal extends ObjectGoal<Database, GSProject> {
	private final ObjectGoal<Database, GSProject> dbGoal;
	private final File dbFile;

	@SafeVarargs
	public LoadDBGoal(GSProject project, ObjectGoal<Database, GSProject> dbGoal, FileGoal<GSProject> updateStoreGoal,
			Goal<GSProject>... dependencies) {
		super(project, GSGoalKey.LOAD_DB,
				project.getDBPath() == null ? append(dependencies, dbGoal, updateStoreGoal) : dependencies);
		this.dbGoal = dbGoal;
		this.dbFile = getProject().getDBPath() == null ? updateStoreGoal.getFile() : new File(getProject().getDBPath());
	}

	@Override
	protected void doMakeThis() {
		try {
			Database db = dbGoal.isMade() ? dbGoal.get()
					: Database.load(dbFile, booleanConfigValue(GSConfigKey.USE_BLOOM_FILTER_FOR_MATCH));
			set(db);
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		} catch (InvalidClassException e) {
			throw new RuntimeException("Database file version does not match genestrip library version.", e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public void setDatabase(Database object) {
		set(object);
	}
}
