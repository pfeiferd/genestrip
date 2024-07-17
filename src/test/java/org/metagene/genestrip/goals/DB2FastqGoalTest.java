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

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.junit.Test;
import org.metagene.genestrip.APITest;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.CountsPerTaxid;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.store.Database;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class DB2FastqGoalTest {
	@Test
	public void testDB2FastqGoal() throws IOException {
		File baseDir = APITest.getBaseDir();

		GSCommon config = new GSCommon(baseDir);

		// Create the 'human_virus' project (whose config files are part of the
		// release).
		GSProject project = new GSProject(config, "human_virus", null, null, null, null, null, false, "64320,12637+",
				null, null, null, false);
		GSMaker maker = new GSMaker(project);

		String[] taxids = new String[] { "64320", "12637", "11053", "11060", "11069", "11070" };
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

		for (int i = 0; i < taxids.length; i++) {
			File file = goal.getOutputFile(taxids[i]);
			MatchingResult result = maker.match(false, "", file.toString());
			Map<String, CountsPerTaxid> map = result.getTaxid2Stats();
			assertEquals(kmers[i], map.get(taxids[i]).getKMers());
			assertEquals(kmers[i], map.get(taxids[i]).getUniqueKMers());
			assertEquals(1, map.size());
		}
		maker.dumpAll();
	}
}
