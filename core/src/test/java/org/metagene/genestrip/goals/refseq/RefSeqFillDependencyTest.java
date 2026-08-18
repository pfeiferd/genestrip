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

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import org.junit.Test;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;

/**
 * Tests that a database built without the RefSeq release does not depend on it either.
 * <p>
 * With {@code refseq.filldb=false} the filling goals read the additional fastas alone - the assemblies
 * fetched from Genbank, say - and never look at a release file. They nevertheless used to carry the
 * download goal among their dependencies, so requesting any of them fetched and checksummed a release
 * of hundreds of gigabytes whose content was then not read. The exception is the goal computing the
 * lowest common ancestors, which reads the release whatever the database was filled from, since the
 * genomes that raise a k-mer above the requested taxa are precisely those of the other taxa.
 */
public class RefSeqFillDependencyTest {
	private static GSProject project(boolean refSeqFillDb) throws IOException {
		File baseDir = new File(System.getProperty("buildDirectory", "target"), "data");
		GSProject project = new GSProject(new GSCommon(baseDir), "human_virus", null, null, null, null, null, null,
				null, null, null, false);
		project.initConfigParam(GSConfigKey.REF_SEQ_DB, refSeqFillDb);
		return project;
	}

	private static boolean dependsOnTheRelease(GSProject project, GSGoalKey key) {
		GSMaker maker = new GSMaker(project);
		@SuppressWarnings("unchecked")
		Goal<GSProject> goal = (Goal<GSProject>) maker.getGoal(key);
		@SuppressWarnings("unchecked")
		Goal<GSProject> download = (Goal<GSProject>) maker.getGoal(GSGoalKey.REFSEQFNA);
		return goal.hasTransDependencyFor(download);
	}

	/** Filling from the release means depending on it, which is the ordinary case. */
	@Test
	public void testTheFillingGoalsDependOnTheReleaseWhenTheyReadIt() throws IOException {
		GSProject project = project(true);
		for (GSGoalKey key : new GSGoalKey[] { GSGoalKey.FILLSIZE, GSGoalKey.TEMPINDEX, GSGoalKey.FILL_DB }) {
			assertTrue(key + " reads the release and must depend on it", dependsOnTheRelease(project, key));
		}
	}

	/** Not filling from it means not fetching it. */
	@Test
	public void testTheFillingGoalsDoNotDependOnTheReleaseWhenTheyIgnoreIt() throws IOException {
		GSProject project = project(false);
		for (GSGoalKey key : new GSGoalKey[] { GSGoalKey.FILLSIZE, GSGoalKey.TEMPINDEX, GSGoalKey.FILL_DB }) {
			assertFalse(key + " does not read the release and must not fetch it", dependsOnTheRelease(project, key));
		}
	}

	/**
	 * The lowest-common-ancestor update is the exception and must keep the dependency either way: a
	 * k-mer that a relative also carries would otherwise stay claimed for the requested taxon.
	 */
	@Test
	public void testTheUpdateDependsOnTheReleaseEitherWay() throws IOException {
		assertTrue(dependsOnTheRelease(project(true), GSGoalKey.UPDATE_DB));
		assertTrue(dependsOnTheRelease(project(false), GSGoalKey.UPDATE_DB));
	}
}
