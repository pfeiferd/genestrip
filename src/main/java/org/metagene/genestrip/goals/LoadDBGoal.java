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

import java.io.*;

import me.tongfei.progressbar.DelegatingProgressBarConsumer;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarBuilder;
import me.tongfei.progressbar.ProgressBarStyle;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;

public class LoadDBGoal<P extends GSProject> extends ObjectGoal<Database, P> implements Goal.LogHeapInfo {
	private final ObjectGoal<Database, P> dbGoal;
	private final File dbFile;

	@SafeVarargs
	public LoadDBGoal(P project, GoalKey key, ObjectGoal<Database, P> dbGoal, FileGoal<P> updateStoreGoal,
					  Goal<P>... dependencies) {
		super(project, key,
				project.getDBPath() == null ? append(dependencies, dbGoal, updateStoreGoal) : dependencies);
		this.dbGoal = dbGoal;
		this.dbFile = getProject().getDBPath() == null ? updateStoreGoal.getFile() : new File(getProject().getDBPath());
	}

	@Override
	protected void doMakeThis() {
		try {
			if (booleanConfigValue(GSConfigKey.PROGRESS_BAR) && !dbGoal.isMade()) {
				try (StreamingResource.StreamAccess sa = new StreamingFileResource(dbFile, true).openStream()) {
					try (ProgressBar pb = GSProgressBarCreator.newGSProgressBar(getKey().getName(), sa, null)) {
						doLoadDB(sa.getInputStream());
					}
				}
			}
			else {
				doLoadDB(StreamProvider.getInputStreamForFile(dbFile));
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected void doLoadDB(InputStream stream) {
		try {
			Database db = dbGoal.isMade() ? dbGoal.get()
					: Database.load(stream, booleanConfigValue(GSConfigKey.USE_BLOOM_FILTER_FOR_MATCH));
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
