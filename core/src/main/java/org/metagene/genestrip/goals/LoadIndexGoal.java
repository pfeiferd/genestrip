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

import me.tongfei.progressbar.ProgressBar;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.bloom.KMerProbFilter;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;

public class LoadIndexGoal<P extends GSProject> extends ObjectGoal<KMerProbFilter, P> implements Goal.LogHeapInfo {
	private final ObjectGoal<KMerProbFilter, P> bloomIndex;
	private final File dbFile;

	@SafeVarargs
	public LoadIndexGoal(P project, ObjectGoal<KMerProbFilter, P> bloomIndex,
			FileGoal<P> storeIndexGoal, Goal<P>... dependencies) {
		super(project, GSGoalKey.LOAD_INDEX,
				project.getDBPath() == null ? append(dependencies, bloomIndex, storeIndexGoal) : dependencies);
		this.bloomIndex = bloomIndex;
		this.dbFile = getProject().getDBPath() == null ? storeIndexGoal.getFile() : new File(getProject().getDBPath());
	}

	@Override
	protected void doMakeThis() {
		try {
			if (booleanConfigValue(GSConfigKey.PROGRESS_BAR) && !bloomIndex.isMade()) {
				try (StreamingResource.StreamAccess sa = new StreamingFileResource(dbFile, false).openStream()) {
					try (ProgressBar pb = GSProgressBarCreator.newGSProgressBar(getKey().getName(), sa, null)) {
						doLoadIndex(sa.getInputStream());
					}
				}
			}
			else {
				doLoadIndex(StreamProvider.getInputStreamForFile(dbFile));
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected void doLoadIndex(InputStream stream) {
		try {
			KMerProbFilter filter = bloomIndex.isMade() ? bloomIndex.get()
					: KMerProbFilter.load(stream);
			set(filter);
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		} catch (InvalidClassException e) {
			throw new RuntimeException("Index file version does not match genestrip library version.", e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
