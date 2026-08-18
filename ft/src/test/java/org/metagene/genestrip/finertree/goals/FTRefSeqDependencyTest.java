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

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import org.junit.Test;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.FinerTreeMaker;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;

/**
 * Tests that the goals building the k-mer index read the RefSeq release whatever the database was
 * filled from.
 * <p>
 * What keeps a k-mer from being pushed down to a child of the node it sits on is a genome that carries
 * it and lies elsewhere. When the release did not fill the database - {@code refseq.filldb=false} -
 * those genomes are precisely the ones only the release holds, so a pass that skipped it would not see
 * them, would set no OTHER slot for them, and would attribute the k-mer to whichever children happen to
 * be in the database: a specificity the genomes outside it contradict.
 */
public class FTRefSeqDependencyTest {
	private static boolean dependsOnTheRelease(boolean refSeqFillDb, GoalKey key) throws IOException {
		File baseDir = new File(System.getProperty("buildDirectory", "target"), "data");
		FTProject project = new FTProject(new GSCommon(baseDir), "human_virus", null, null, null, null, null, null,
				null, null, null, false);
		project.initConfigParam(GSConfigKey.REF_SEQ_DB, refSeqFillDb);
		FinerTreeMaker<FTProject> maker = new FinerTreeMaker<>(project);
		@SuppressWarnings("unchecked")
		Goal<FTProject> goal = (Goal<FTProject>) maker.getGoal(key);
		@SuppressWarnings("unchecked")
		Goal<FTProject> download = (Goal<FTProject>) maker.getGoal(GSGoalKey.REFSEQFNA);
		return goal.hasTransDependencyFor(download);
	}

	/** The filter that the reassignment reads must be built from every genome, not only the database's. */
	@Test
	public void testTheIndexGoalsDependOnTheReleaseEitherWay() throws IOException {
		for (GoalKey key : new GoalKey[] { FTGoalKey.KMER_INDEX_BLOOM, FTGoalKey.KMER_INDEX_SIZE }) {
			assertTrue(key + " must read the release when the database was filled from it",
					dependsOnTheRelease(true, key));
			assertTrue(key + " must read the release even when the database was not filled from it",
					dependsOnTheRelease(false, key));
		}
	}
}
