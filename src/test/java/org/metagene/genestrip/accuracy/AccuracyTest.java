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
import org.metagene.genestrip.GSConfig;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileDownloadGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.tax.TaxTree;

@Ignore
public class AccuracyTest {
	private final GSConfig config;
	private final GSProject project;
	private final GSMaker maker;
	private final FastaTransformGoal fastaTransformGoal;
	private final AccuracyProjectGoal accuracyProjectGoal;
	private final TaxIdsTxtGoal taxIdsTxtGoal;
	private final File multiMatchCSVFile;

	@SuppressWarnings("unchecked")
	public AccuracyTest() {
		try {
			config = new GSConfig(getTargetDir()) {
				// Override for fair speed test.
				@Override
				public int getThreads() {
					return 0;
				}
			};
			project = new GSProject(config, "accuracy", 31, null, null, null, null);
			maker = new GSMaker(project);

			multiMatchCSVFile = new File(project.getProjectDir(), "multimatch.txt");

			AccuracyDataDownloadGoal accDownloadGoal = new AccuracyDataDownloadGoal(project, "accdatadownload",
					maker.getGoal("setup"));

			FileDownloadGoal<GSProject> acc2taxidDownloadGoal = new FileDownloadGoal<GSProject>(project, "accdownload",
					maker.getGoal("commonsetup")) {
				@Override
				public List<File> getFiles() {
					return Arrays.asList(
							new File(project.getConfig().getCommonDir(), AccessionNumber2TaxidGoal.ACCESSION_MAP_FILE),
							new File(project.getConfig().getCommonDir(),
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
					project, "acc2taxid", new File(project.getFastaDir(), "accuracy/simBA5_accuracy.fa"),
					acc2taxidDownloadGoal);

			accuracyProjectGoal = new AccuracyProjectGoal(project, "accproject");
			taxIdsTxtGoal = new TaxIdsTxtGoal(project, "taxidtxt",
					(ObjectGoal<TaxTree, GSProject>) maker.getGoal("taxtree"), accessionNumber2TaxidGoal);

			fastaTransformGoal = new FastaTransformGoal(project, "fastatransform", accessionNumber2TaxidGoal,
					accDownloadGoal);
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

		AccuracyMatchGoal tempdbMatchGoal = new AccuracyMatchGoal(project, "acctempdbmatch", multiMatchCSVFile,
				(ObjectGoal<TaxTree, GSProject>) maker.getGoal("taxtree"),
				(ObjectGoal<KMerStoreWrapper, GSProject>) maker.getGoal("filleddb"), fastaTransformGoal);
		tempdbMatchGoal.make();

		AccuracyMatchGoal matchGoal = new AccuracyMatchGoal(project, "accmatch", multiMatchCSVFile,
				(ObjectGoal<TaxTree, GSProject>) maker.getGoal("taxtree"),
				(ObjectGoal<KMerStoreWrapper, GSProject>) maker.getGoal("updateddb"), fastaTransformGoal);
		matchGoal.make();
	}

	@Test
	@Ignore
	@SuppressWarnings("unchecked")
	public void testKrakenAccuracy() throws IOException {
		accuracyProjectGoal.make();
		taxIdsTxtGoal.make();
		fastaTransformGoal.make();

		KrakenAccuracyMatchGoal krakenAccuracyMatchGoal = new KrakenAccuracyMatchGoal(project, "krakenacc",
				multiMatchCSVFile, (ObjectGoal<TaxTree, GSProject>) maker.getGoal("taxtree"), fastaTransformGoal);
		krakenAccuracyMatchGoal.make();
	}
}
