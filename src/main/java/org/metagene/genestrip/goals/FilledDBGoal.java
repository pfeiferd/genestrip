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

import java.io.IOException;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;

public class FilledDBGoal extends ObjectGoal<Database, GSProject> {
	private final ObjectGoal<Database, GSProject> fillDBGoal;
	private final FileGoal<GSProject> filledStoreGoal;

	@SafeVarargs
	public FilledDBGoal(GSProject project, ObjectGoal<Database, GSProject> fillDBGoal,
			FileGoal<GSProject> filledStoreGoal, Goal<GSProject>... dependencies) {
		super(project, GSGoalKey.LOAD_TEMPDB, append(dependencies, filledStoreGoal, fillDBGoal));
		this.fillDBGoal = fillDBGoal;
		this.filledStoreGoal = filledStoreGoal;
	}

	@Override
	protected void doMakeThis() {
		try {
			Database db = fillDBGoal.isMade() ? fillDBGoal.get() : Database.load(filledStoreGoal.getFile(), true);
			set(db);
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
