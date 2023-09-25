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
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AccessionFileProcessor;
import org.metagene.genestrip.refseq.RefSeqCategory;

public class AccessionMapSizeGoal extends ObjectGoal<Integer, GSProject> {
	private final ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal;
	private final RefSeqCatalogDownloadGoal catalogGoal;

	@SafeVarargs
	public AccessionMapSizeGoal(GSProject project, String name,
			ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal, RefSeqCatalogDownloadGoal catalogGoal,
			RefSeqFnaFilesDownloadGoal downloadGoal, Goal<GSProject>... deps) {
		super(project, name, Goal.append(deps, categoriesGoal, catalogGoal, downloadGoal));
		this.categoriesGoal = categoriesGoal;
		this.catalogGoal = catalogGoal;
	}

	@Override
	public void makeThis() {
		File catalogFile = catalogGoal.getCatalogFile();
		AccessionFileProcessor processor = new AccessionFileProcessor(categoriesGoal.get()[1],
				getProject().isUseCompleteGenomesOnly()) {
			private int counter = 0;

			@Override
			public void processCatalog(File catalogFile) {
				super.processCatalog(catalogFile);
				set(counter);

				if (getLogger().isInfoEnabled()) {
					getLogger().info("Map size determined: " + counter);
				}
			}

			@Override
			protected void handleEntry(byte[] target, int taxIdEnd, int accessionStart, int accessionEnd) {
				counter++;
			}
		};
		processor.processCatalog(catalogFile);
	}
}
