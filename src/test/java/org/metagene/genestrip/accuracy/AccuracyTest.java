package org.metagene.genestrip.accuracy;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.junit.Test;
import org.metagene.genestrip.GSConfig;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileDownloadGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.tax.TaxTree;

// @Ignore
public class AccuracyTest {
	@SuppressWarnings("unchecked")
	@Test
	public void testAccuracy() throws IOException {
		File testBaseDir = getTargetDir();
		GSConfig config = new GSConfig(testBaseDir);
		GSProject project = new GSProject(config, "accuracy", 31, null, null, null, null);

		GSMaker maker = new GSMaker(project);

		AccuracyDataDownloadGoal accDownloadGoal = new AccuracyDataDownloadGoal(project, "accdatadownload",
				maker.getGoal("setup"));

		FileDownloadGoal<GSProject> acc2taxidDownloadGoal = new FileDownloadGoal<GSProject>(project, "accdownload",
				maker.getGoal("commonsetup")) {
			@Override
			public List<File> getFiles() {
				return Arrays.asList(
						new File(project.getConfig().getCommonDir(), AccessionNumber2TaxidGoal.ACCESSION_MAP_FILE),
						new File(project.getConfig().getCommonDir(), AccessionNumber2TaxidGoal.OLD_ACCESSION_MAP_FILE));
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

		ObjectGoal<Map<String, String>, GSProject> accessionNumber2TaxidGoal = new AccessionNumber2TaxidGoal(project,
				"acc2taxid", new File(project.getFastaDir(), "accuracy/simBA5_accuracy.fa"), acc2taxidDownloadGoal);

		new TaxIdsTxtGoal(project, "taxidtxt", (ObjectGoal<TaxTree, GSProject>) maker.getGoal("taxtree"),
				accessionNumber2TaxidGoal, new AccuracyProjectGoal(project, "accproject")).make();

		FastaTransformGoal fastaTransformGoal = new FastaTransformGoal(project, "fastatransform",
				accessionNumber2TaxidGoal, accDownloadGoal);

		// fastaTransformGoal.make();

		AccuracyMatchGoal matchGoal = new AccuracyMatchGoal(project, "accmatch",
				new File(project.getProjectDir(), "multimatch.txt"),
				(ObjectGoal<TaxTree, GSProject>) maker.getGoal("taxtree"),
				(ObjectGoal<KMerStoreWrapper, GSProject>) maker.getGoal("updateddb"), fastaTransformGoal);

		matchGoal.make();
	}

	protected File getTargetDir() {
		String projectDir = System.getProperty("project.build.directory");
		if (projectDir != null) {
			return new File(projectDir);
		}
		String relPath = getClass().getProtectionDomain().getCodeSource().getLocation().getFile();
		return new File(relPath).getParentFile();
	}
}
