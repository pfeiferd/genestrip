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
package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerStoreWrapper;

public class FilledStoreGoal extends ObjectGoal<KMerStoreWrapper, GSProject> {
	private final FileGoal<GSProject> fillStoreGoal;

	@SafeVarargs
	public FilledStoreGoal(GSProject project, String name, FileGoal<GSProject> fillStoreGoal,
			Goal<GSProject>... dependencies) {
		super(project, name, Goal.append(dependencies, fillStoreGoal));
		this.fillStoreGoal = fillStoreGoal;
	}

	@Override
	public void makeThis() {
		try {
			File filterFile = null;
			if (getProject().isUseBloomFilterForMatch()) {
				filterFile = getProject().getFilterFile(fillStoreGoal);
				if (!filterFile.exists()) {
					if (getLogger().isWarnEnabled()) {
						getLogger().warn("Missing filter file " + filterFile + ". The database will be used without it.");
					}
					filterFile = null;
				}
			}
			KMerStoreWrapper wrapper = KMerStoreWrapper.load(fillStoreGoal.getFile(), filterFile);
			set(wrapper);
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public void setStoreWrapper(KMerStoreWrapper object) {
		set(object);
	}
}
