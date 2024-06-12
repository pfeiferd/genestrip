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
package org.metagene.genestrip.accuracy;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.junit.Ignore;
import org.junit.Test;
import org.metagene.genestrip.DefaultExecutionContext;
import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.goals.GSFileDownloadGoal;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.make.FileDownloadGoal;
import org.metagene.genestrip.make.GoalKey.DefaultGoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.TaxTree;

@Ignore
public class AccuracyTest {
	private final GSCommon config;
	private final GSProject project;
	private final GSMaker maker;
	private final FastaTransformGoal fastaTransformGoal;
	private final AccuracyProjectGoal accuracyProjectGoal;
	private final TaxIdsTxtGoal taxIdsTxtGoal;

	@SuppressWarnings("unchecked")
	public AccuracyTest() {
		try {
			config = new GSCommon(getTargetDir());
			project = new GSProject(config, "accuracy", "multimatch.txt");
			// Override for fair speed test.
			project.initConfigParam(GSConfigKey.THREADS, 0);
			maker = new GSMaker(project);

			AccuracyDataDownloadGoal accDownloadGoal = new AccuracyDataDownloadGoal(project,
					maker.getGoal(GSGoalKey.SETUP));

			FileDownloadGoal<GSProject> acc2taxidDownloadGoal = new GSFileDownloadGoal(project,
					new DefaultGoalKey("acc2taxiddownload"), maker.getGoal(GSGoalKey.COMMON_SETUP)) {
				@Override
				public List<File> getFiles() {
					return Arrays.asList(
							new File(project.getCommon().getCommonDir(), AccessionNumber2TaxidGoal.ACCESSION_MAP_FILE),
							new File(project.getCommon().getCommonDir(),
									AccessionNumber2TaxidGoal.OLD_ACCESSION_MAP_FILE));
				}

				@Override
				protected boolean isUseHttp() {
					return true;
				}

				protected String getHttpBaseURL() {
					return "https://ftp.ncbi.nih.gov";
				}

				@Override
				protected String getFTPDir(File file) {
					return "/pub/taxonomy/accession2taxid";
				}
			};

			ObjectGoal<Map<String, String>, GSProject> accessionNumber2TaxidGoal = new AccessionNumber2TaxidGoal(
					project, new File(project.getFastaDir(), "accuracy/simBA5_accuracy.fa"),
					acc2taxidDownloadGoal);

			accuracyProjectGoal = new AccuracyProjectGoal(project);
			taxIdsTxtGoal = new TaxIdsTxtGoal(project,
					(ObjectGoal<TaxTree, GSProject>) maker.getGoal(GSGoalKey.TAXTREE), accessionNumber2TaxidGoal);

			fastaTransformGoal = new FastaTransformGoal(project, accessionNumber2TaxidGoal, accDownloadGoal);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected File getTargetDir() {
		String projectDir = System.getProperty("project.build.directory");
		if (projectDir != null) {
			return new File(projectDir);
		}
		String relPath = getClass().getProtectionDomain().getCodeSource().getLocation().getFile();
		return new File(relPath).getParentFile();
	}

	@Test
	@SuppressWarnings("unchecked")
	public void testGenestripAccuracy() throws IOException {
		accuracyProjectGoal.make();
		taxIdsTxtGoal.make();
		fastaTransformGoal.make();

		ExecutionContext bundle = new DefaultExecutionContext(0,
				project.longConfigValue(GSConfigKey.LOG_PROGRESS_UPDATE_CYCLE));

		ObjectGoal<Map<String, List<StreamingResource>>, GSProject> fastqMapGoal = (ObjectGoal<Map<String, List<StreamingResource>>, GSProject>) maker
				.getGoal(GSGoalKey.FASTQ_MAP);

		ObjectGoal<TaxTree, GSProject> taxTreeGoal = (ObjectGoal<TaxTree, GSProject>) maker.getGoal(GSGoalKey.TAXTREE);

		AccuracyMatchGoal tempdbMatchGoal = new AccuracyMatchGoal(project, new DefaultGoalKey("acctempdbmatch"),
				fastqMapGoal, (ObjectGoal<Database, GSProject>) maker.getGoal(GSGoalKey.LOAD_TEMPDB), taxTreeGoal, bundle,
				fastaTransformGoal);
		tempdbMatchGoal.make();

		AccuracyMatchGoal matchGoal = new AccuracyMatchGoal(project, new DefaultGoalKey("accmatch"), fastqMapGoal,
				(ObjectGoal<Database, GSProject>) maker.getGoal(GSGoalKey.LOAD_DB), taxTreeGoal, bundle,
				fastaTransformGoal);
		matchGoal.make();
	}

	@Test
	@Ignore
	@SuppressWarnings("unchecked")
	public void testKrakenAccuracy() throws IOException {
		accuracyProjectGoal.make();
		taxIdsTxtGoal.make();
		fastaTransformGoal.make();

		ObjectGoal<Map<String, List<StreamingResource>>, GSProject> fastqMapGoal = (ObjectGoal<Map<String, List<StreamingResource>>, GSProject>) maker
				.getGoal(GSGoalKey.FASTQ_MAP);

		KrakenAccuracyMatchGoal krakenAccuracyMatchGoal = new KrakenAccuracyMatchGoal(project,
				fastqMapGoal, (ObjectGoal<TaxTree, GSProject>) maker.getGoal(GSGoalKey.TAXTREE), fastaTransformGoal);
		krakenAccuracyMatchGoal.make();
	}
}
