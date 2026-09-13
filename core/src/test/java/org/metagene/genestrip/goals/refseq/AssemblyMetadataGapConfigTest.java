/*
 * Genestrip
 */
package org.metagene.genestrip.goals.refseq;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.File;

import org.junit.Test;
import org.metagene.genestrip.APITest;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;

/**
 * Tests that the assembly metadata gap really is a configuration value.
 * <p>
 * {@link AssemblySizeIndexTest} covers what a gap <i>means</i> - which assembly a sum of sequence
 * lengths is allowed to match. It does so by passing the gap to
 * {@link org.metagene.genestrip.refseq.AssemblySizeIndex#nearest(int, long, int)} directly, so it
 * would still pass if the goal ignored the configuration and used a literal. These tests close that
 * hole: they check the default, that an override is picked up, and that the key keeps its name.
 */
public class AssemblyMetadataGapConfigTest {

	/** Returns the assembly metadata goal of a freshly made project. */
	private AssemblyMetadataGoal<GSProject> goalOf(GSProject project) {
		GSMaker<GSProject> maker = new GSMaker<>(project);
		return (AssemblyMetadataGoal<GSProject>) maker.getGoal(GSGoalKey.ASSEMBLYMETA);
	}

	private GSProject project() throws Exception {
		File baseDir = APITest.getBaseDir();
		return new GSProject(new GSCommon(baseDir), "human_virus", true);
	}

	/**
	 * The documented default is 100 bases. It is not a free choice: measurements put the share of
	 * multi-replicon assemblies reconstructed at 74.8% for gap 0, 82.1% for gap 10 and 92.0% for gap
	 * 100, so lowering it silently would cost recall.
	 */
	@Test
	public void testDefaultGapIs100() throws Exception {
		assertEquals(100, goalOf(project()).getGap());
	}

	/** An override in the project configuration must reach the goal. */
	@Test
	public void testGapIsConfigurable() throws Exception {
		GSProject project = project();
		project.initConfigParam(GSConfigKey.ASSEMBLY_METADATA_GAP, 250);
		assertEquals(250, goalOf(project).getGap());
	}

	/**
	 * The goal must read the value when it runs rather than latch it at construction, so that a
	 * configuration applied after the goal graph is built still takes effect.
	 */
	@Test
	public void testGapIsNotLatchedAtConstruction() throws Exception {
		GSProject project = project();
		AssemblyMetadataGoal<GSProject> goal = goalOf(project);
		assertEquals(100, goal.getGap());
		project.initConfigParam(GSConfigKey.ASSEMBLY_METADATA_GAP, 7);
		assertEquals("the goal latched the gap instead of reading it", 7, goal.getGap());
	}

	/**
	 * The name is what a user writes in a project's config.properties, so a rename is a breaking
	 * change and should have to be made here as well.
	 */
	@Test
	public void testGapKeyName() {
		assertEquals("refseq.assemblyMetadataGap", GSConfigKey.ASSEMBLY_METADATA_GAP.getName());
	}

	/**
	 * A negative gap is meaningless and the key's declared range rejects it: the assignment fails and
	 * the previous value stands, rather than a negative gap reaching the index and matching nothing.
	 */
	@Test
	public void testNegativeGapIsRejected() throws Exception {
		GSProject project = project();
		assertFalse("a negative gap must not be accepted",
				project.initConfigParam(GSConfigKey.ASSEMBLY_METADATA_GAP, -1));
		assertEquals("a rejected gap must leave the default standing", 100, goalOf(project).getGap());
	}

	/** Gap 0 is a legitimate setting - it means the sum must match an assembly's length exactly. */
	@Test
	public void testGapZeroIsAllowed() throws Exception {
		GSProject project = project();
		assertTrue(project.initConfigParam(GSConfigKey.ASSEMBLY_METADATA_GAP, 0));
		assertEquals(0, goalOf(project).getGap());
	}

	/** The goal is registered under its key and is of the expected type. */
	@Test
	public void testGoalIsRegistered() throws Exception {
		assertNotNull(goalOf(project()));
	}
}
