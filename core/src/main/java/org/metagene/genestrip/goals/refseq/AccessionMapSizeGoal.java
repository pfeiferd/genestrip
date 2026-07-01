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
import java.util.List;
import java.util.Set;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSConfigKey.SeqType;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AccessionFileProcessor;
import org.metagene.genestrip.refseq.RefSeqCategory;

/**
 * Goal that counts how many accession entries the RefSeq catalog contributes for the selected
 * categories and sequence type; the resulting size is used to pre-allocate the
 * {@link org.metagene.genestrip.refseq.AccessionMap}.
 *
 * @param <P> the project type
 */
public class AccessionMapSizeGoal<P extends GSProject> extends ObjectGoal<Integer, P> {
	private final ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal;
	private final RefSeqCatalogDownloadGoal catalogGoal;

	/**
	 * Creates the goal, wiring the categories and catalog download goals it scans.
	 *
	 * @param project        the project type
	 * @param categoriesGoal the goal providing the selected RefSeq categories
	 * @param catalogGoal    the goal providing the downloaded RefSeq catalog file
	 * @param deps           additional goals this goal depends on
	 */
	@SafeVarargs
	public AccessionMapSizeGoal(P project,
			ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal, RefSeqCatalogDownloadGoal catalogGoal,
			Goal<P>... deps) {
		super(project, GSGoalKey.ACCMAPSIZE, Goal.append(deps, categoriesGoal, catalogGoal));
		this.categoriesGoal = categoriesGoal;
		this.catalogGoal = catalogGoal;
	}

	@Override
	protected void doMakeThis() {
		File catalogFile = catalogGoal.getCatalogFile();
		AccessionFileProcessor processor = new AccessionFileProcessor(categoriesGoal.get(),
				(SeqType) configValue(GSConfigKey.SEQ_TYPE),
				(List<GSConfigKey.RefSeqStatus>) configValue(GSConfigKey.RES_SEQ_STATUS)) {
			private int counter = 0;

			@Override
			public void processCatalog(StreamingResource catalogFile) {
				super.processCatalog(catalogFile);
				set(counter);

				if (getLogger().isDebugEnabled()) {
					getLogger().debug("Map size determined: " + counter);
				}
			}

			@Override
			protected boolean isProgressBar() {
				return booleanConfigValue(GSConfigKey.PROGRESS_BAR);
			}

			protected String getProgressBarTaskName() {
				return getKey().getName();
			}

			@Override
			protected void handleEntry(byte[] target, int taxIdEnd, int accessionStart, int accessionEnd) {
				counter++;
			}
		};
		processor.processCatalog(new StreamingFileResource(catalogFile));
	}
}
