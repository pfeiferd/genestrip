package org.metagene.genestrip.goals;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.junit.Test;
import org.metagene.genestrip.APITest;
import org.metagene.genestrip.GSConfig;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.CountsPerTaxid;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.util.DefaultExecutorServiceBundle;
import org.metagene.genestrip.util.ExecutorServiceBundle;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class DB2FastqGoalTest {
	@Test
	public void testDB2FastqGoal() throws IOException {
		File baseDir = APITest.getBaseDir();

		GSConfig config = new GSConfig(baseDir);

		// Create the 'human_virus' project (whose config files are part of the
		// release).
		GSProject project = new GSProject(config, "human_virus", 31, null, null, null, null, "64320,12637+");
		GSMaker maker = new GSMaker(project);

		String[] taxids = new String[] { "64320", "12637", "11053", "11060", "11069", "11070" };
		@SuppressWarnings("unchecked")
		ObjectGoal<KMerStoreWrapper, GSProject> storeGoal = (ObjectGoal<KMerStoreWrapper, GSProject>) maker
				.getGoal("updateddb");
		Object2LongMap<String> stats = storeGoal.get().getStats();
		long[] kmers = new long[taxids.length];
		for (int i = 0; i < taxids.length; i++) {
			kmers[i] = stats.getLong(taxids[i]);
		}

		DB2FastqGoal goal = (DB2FastqGoal) maker.getGoal("db2fastq");
		goal.cleanThis();
		goal.make();

		int[] count = new int[1];
		for (int i = 0; i < taxids.length; i++) {
			count[0] = i;
			File file = goal.getOutputFile(taxids[i]);
			project = new GSProject(config, "human_virus", 31, null, file, null, null, null);
			maker = new GSMaker(project);

			ExecutorServiceBundle bundle = new DefaultExecutorServiceBundle(0, config.getMatchLogUpdateCycle());
			@SuppressWarnings("unchecked")
			Goal<GSProject> myMatchGoal = new MultiMatchGoal(project, "mymatch", false, file,
					(ObjectGoal<TaxTree, GSProject>) maker.getGoal("taxtree"),
					(ObjectGoal<KMerStoreWrapper, GSProject>) maker.getGoal("updateddb"), false, bundle) {
				protected void writeOutputFile(File file, MatchingResult result, KMerStoreWrapper wrapper)
						throws IOException {
					Map<String, CountsPerTaxid> map = result.getTaxid2Stats();
					int i = count[0];
					assertEquals(kmers[i], map.get(taxids[i]).getKMers());
					assertEquals(kmers[i], map.get(taxids[i]).getUniqueKMers());
					assertEquals(1, map.size());
				}
			};
			myMatchGoal.make();
		}
	}
}
