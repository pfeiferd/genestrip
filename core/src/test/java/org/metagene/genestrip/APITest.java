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
package org.metagene.genestrip;

import java.io.File;
import java.io.IOException;

import org.junit.BeforeClass;
import org.junit.Test;
import org.metagene.genestrip.make.Goal;

/**
 * Demonstrates (and tests) how to drive Genestrip programmatically through its API, using the
 * bundled {@code human_virus} sample project.
 * <p>
 * Two equivalent styles are shown:
 * <ul>
 * <li>{@link #testAndDemoSimpleAPI()} uses the high-level convenience methods
 *     {@link GSMaker#match} / {@link GSMaker#filter}, while</li>
 * <li>{@link #testAndDemoBasicAPI()} uses the lower-level goal API directly
 *     ({@link GSMaker#getGoal}, {@link Goal#cleanThis()}, {@link Goal#make()}).</li>
 * </ul>
 * In both cases, making a {@code match} or {@code filter} goal automatically pulls in whatever
 * upstream goals are still missing - e.g. the {@code db} goal that builds the database - so the
 * examples run end to end even on a fresh checkout.
 */
public class APITest {
	// Runs once before all tests in this class: run the 'clear' goal, which empties the project's
	// generated-output folders (csv, log and krakenout) so each test starts from a clean slate.
	// Note: 'clear' does not delete the database itself; the match/filter goals below rebuild their
	// own results and (re)create the database only if it is missing.
	@BeforeClass
	public static void clearDB() throws IOException {
		File baseDir = getBaseDir();

		// Load Genestrip's system-wide configuration, rooted at the base directory.
		GSCommon config = new GSCommon(baseDir);

		// Create the 'human_virus' project (whose config files ship as part of the release, under
		// the base directory).
		GSProject project = new GSProject(config, "human_virus");
		GSMaker<GSProject> maker = new GSMaker<>(project);
		maker.getGoal(GSGoalKey.CLEAR).make();

		// Release the resources held by the maker (data cached in memory and its worker threads).
		maker.dumpAll();
	}

	// Demonstrates the high-level convenience API: GSMaker.match(...) / GSMaker.filter(...).
	@Test
	public void testAndDemoSimpleAPI() throws IOException {
		File baseDir = getBaseDir();

		// Load Genestrip's system-wide configuration, rooted at the base directory.
		GSCommon config = new GSCommon(baseDir);

		// Create the 'human_virus' project (whose config files ship as part of the release).
		GSProject project = new GSProject(config, "human_virus");
		GSMaker<GSProject> maker = new GSMaker<>(project);

		// The fastq file to run a match on (here a small sample file shipped with the release).
		File fastq = new File(getBaseDir(), "projects/human_virus/fastq/sample.fastq.gz");
		// Example of setting a configuration parameter programmatically (this overrides the value
		// from the project's config files); see GSConfigKey for the available keys.
		project.initConfigParam(GSConfigKey.USE_BLOOM_FILTER_FOR_MATCH, true);
		//project.initConfigParam(GSConfigKey.MAX_KMER_RES_COUNTS, 35);
		// Run the 'match' goal on the fastq file. Arguments:
		//   1st boolean: use the 'matchlr' variant (matching without read classification, intended
		//                for long reads) when true; 'false' selects the standard 'match'.
		//   2nd boolean: when true, first delete any previously generated result for this input.
		//   key (null): optional name for the input group; null uses a default derived from the file.
		// Making 'match' also triggers any prerequisite goals, e.g. the 'db' goal that builds the
		// 'human_virus' database if it does not exist yet.
		maker.match(false, true, null, fastq.toString());

		// Run the 'filter' goal on the same input. The leading 'true' first deletes any previously
		// generated result. Making 'filter' triggers any prerequisite goals, e.g. the 'index' goal
		// that builds the 'human_virus' filtering database if needed.
		maker.filter(true, null, fastq.toString());

		// Release the resources held by the maker (data cached in memory and its worker threads).
		maker.dumpAll();
	}

	// Demonstrates the lower-level goal API: obtain a Goal and drive it with cleanThis()/make().
	@Test
	public void testAndDemoBasicAPI() throws IOException {
		// Get the base directory for Genestrip. (Depends on your setup.)
		File baseDir = getBaseDir();

		// Load Genestrip's system-wide configuration, rooted at the base directory.
		GSCommon config = new GSCommon(baseDir);

		// The fastq file to run a match on (here a small sample file shipped with the release).
		File fastq = new File(getBaseDir(), "projects/human_virus/fastq/sample.fastq.gz");

		// Create the 'human_virus' project and bind the fastq input to it up front (the last
		// constructor argument). The goals below then read their input from the project, so - unlike
		// the simple API above - no file paths need to be passed per call. The 'key' argument (here
		// null) optionally names the input group.
		GSProject project = new GSProject(config, "human_virus", null, new String[] { fastq.getCanonicalPath() });
		GSMaker<GSProject> maker = new GSMaker<>(project);
		// Obtain the 'match' goal. Making it triggers any prerequisite goals, e.g. the 'db' goal that
		// builds the 'human_virus' database if it does not exist yet.
		Goal<GSProject> goal = maker.getGoal(GSGoalKey.MATCH);
		// First, delete a potential result that has previously been generated.
		goal.cleanThis();
		// then run the match.
		goal.make();

		// Obtain the 'filter' goal. Making it triggers any prerequisite goals, e.g. the 'index' goal
		// that builds the 'human_virus' filtering database if needed.
		goal = maker.getGoal(GSGoalKey.FILTER);
		// First, delete a potential result that has previously been generated.
		goal.cleanThis();
		// then run the filter.
		goal.make();

		// Release the resources held by the maker (data cached in memory and its worker threads).
		maker.dumpAll();
	}

	// Locates the base directory of the bundled sample projects that ship as part of the release.
	// (For your own projects this is typically a fixed path you choose instead.)
	public static File getBaseDir() {
		return new File(getProjectDir(), "data");
	}

	public static File getProjectDir() {
		String projectDir = System.getProperty("maven.multiModuleProjectDirectory");
		if (projectDir != null) {
			return new File(projectDir);
		}
		String relPath = APITest.class.getProtectionDomain().getCodeSource().getLocation().getFile();
		return new File(relPath).getParentFile().getParentFile().getParentFile();
	}
}
