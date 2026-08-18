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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.nio.file.Files;
import java.util.List;
import java.util.Arrays;
import java.util.Map;

import org.junit.Test;
import org.metagene.genestrip.APITest;
import org.metagene.genestrip.goals.refseq.DBGoalTest;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey.DefaultGoalKey;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.goals.refseq.ComprehensiveFilterTest;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.CountsPerTaxid;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.store.Database;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class DB2FastqGoalTest extends ComprehensiveFilterTest {
	@Override
	public void testUpdate() throws IOException {
		// Just to avoid running the test from the superclass ...
	}

	@Override
	public void testKrakenOutput() throws IOException {
		// Just to avoid running the test from the superclass ...
	}

	@Override
	public <P extends GSProject> void testFilterOutput() throws IOException {
		// Just to avoid running the test from the superclass ...
	}

	@Test
	public void testDB2FastqGoal() throws IOException {
		// The project gets a copy of its own below the build output, so that this test neither reads
		// nor writes what the other tests of the release's `human_virus' do - they configure it
		// differently, and whichever ran last would decide what this one finds in its database.
		// Everything that is downloaded rather than configured stays shared: the common folder is the
		// release's, since re-fetching the RefSeq catalogue and the viral genomes for a test would
		// cost gigabytes.
		File releaseCommon = new File(APITest.getBaseDir(), "common");
		GSCommon config = new GSCommon(getBaseDir()) {
			@Override
			public File getCommonDir() {
				return releaseCommon;
			}
		};

		String[] taxids = new String[] { "64320", "12637", "11053", "11060", "11069", "11070" };

		// Create the 'human_virus' project, whose configuration files are copied in below.
		GSProject project = new GSProject(config, "human_virus", null, null, null, null, null, "64320,12637+",
				null, null, null, false);
		project.initConfigParam(GSConfigKey.GZIP_FASTQ_OUTPUT, false);
		project.initConfigParam(GSConfigKey.TAX_IDS, Arrays.asList(taxids));

		new HumanVirusProjectGoal(project).make();

		GSMaker maker = new GSMaker(project);

		@SuppressWarnings("unchecked")
		ObjectGoal<Database, GSProject> storeGoal = (ObjectGoal<Database, GSProject>) maker.getGoal(GSGoalKey.LOAD_DB);
		Object2LongMap<String> stats = storeGoal.get().getStats();
		long[] kmers = new long[taxids.length];
		for (int i = 0; i < taxids.length; i++) {
			kmers[i] = stats.getLong(taxids[i]);
		}

		DB2FastqGoal goal = (DB2FastqGoal) maker.getGoal(GSGoalKey.DB2FASTQ);
		goal.cleanThis();
		goal.make();

		// int[] mapSizes = new int[] { 12, 12, 13, 13, 13, 13 };
 		for (int i = 0; i < taxids.length; i++) {
			File file = goal.getOutputFile(taxids[i]);
			MatchingResult result = maker.cleanMatch(false, taxids[i], file.toString());
			Map<String, CountsPerTaxid> map = result.getTaxid2Stats();

			assertEquals(kmers[i], map.get(taxids[i]).getKMers());
			assertEquals(kmers[i], map.get(taxids[i]).getUniqueKMers());
			// System.out.println(map.size());
			// Depends on RefSeq version, apparently:
			System.out.println(map.size());
			assertTrue(12 == map.size() || 13 == map.size() || 14 == map.size());
		}
		maker.dumpAll();
	}

	/**
	 * Puts the configuration files of the {@code human_virus} project in place below the build output,
	 * copying them from the test resources exactly as {@link DBGoalTest.Dengue1ProjectGoal} does for
	 * its own project. Only the files that describe the project are copied; everything else is derived
	 * from them.
	 */
	public static class HumanVirusProjectGoal extends FileListGoal<GSProject> {
		@SafeVarargs
		public HumanVirusProjectGoal(GSProject project, Goal<GSProject>... dependencies) {
			super(project, new DefaultGoalKey("human_virus_project"), (List<File>) null, dependencies);
			addFile(new File(project.getProjectDir(), "taxids.txt"));
			addFile(new File(project.getProjectDir(), "categories.txt"));
		}

		@Override
		protected void makeFile(File file) throws IOException {
			// Every directory the maker's own setup goal would create, not just the ones holding the
			// configuration: that goal creates them with mkdir(), which fails silently when its parent
			// is missing, and a clean-and-remake cycle within one maker does not run it a second time.
			for (File dir : new File[] { getProject().getCommon().getBaseDir(), getProject().getProjectsDir(),
					getProject().getProjectDir(), getProject().getFastaDir(), getProject().getFastqDir(),
					getProject().getDBDir(), getProject().getKrakenOutDir(), getProject().getResultsDir(),
					getProject().getLogDir(), getProject().getFastqResDir() }) {
				if (!dir.exists()) {
					dir.mkdirs();
				}
			}
			if (!file.exists()) {
				URL resource = getClass().getClassLoader().getResource("projects/human_virus/" + file.getName());
				if (resource == null) {
					throw new IOException("Missing test resource projects/human_virus/" + file.getName());
				}
				Files.copy(new File(resource.getFile()).toPath(), file.toPath());
			}
		}
	}
}
