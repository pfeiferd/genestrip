package org.metagene.genestrip.accuracy;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;
import org.metagene.genestrip.GSConfig;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileDownloadGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.tax.TaxTree;

public class AccuracyTest {
	@Test
	public void testAccuracy() throws IOException {
		File testBaseDir = getTargetDir();
		GSConfig config = new GSConfig(testBaseDir);
		GSProject project = new GSProject(config, "accuracy", 31, null, null, null, null);

		GSMaker maker = new GSMaker(project);

		new AccuracyProjectGoal(project, "accproject").make();

		AccuracyDataDownloadGoal accDownloadGoal = new AccuracyDataDownloadGoal(project, "accdatadownload",
				maker.getGoal("setup"));
		
		FileDownloadGoal<GSProject> acc2taxidDownloadGoal = new FileDownloadGoal<GSProject>(project, "assdownload", maker.getGoal("commonsetup")) {
			@Override
			public List<File> getFiles() {
				return Arrays.asList(new File(project.getConfig().getCommonDir(), FastaTransformGoal.ACCESSION_MAP_FILE),
						new File(project.getConfig().getCommonDir(), FastaTransformGoal.OLD_ACCESSION_MAP_FILE));
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
		
		FastaTransformGoal fastaTransformGoal = new FastaTransformGoal(project, "fastatransform", accDownloadGoal, acc2taxidDownloadGoal);
		
		fastaTransformGoal.make();

		@SuppressWarnings("unchecked")
		AccuracyMatchGoal matchGoal = new AccuracyMatchGoal(project, "accmatch",
				new File(project.getProjectDir(), "multimatch.txt"),
				(ObjectGoal<TaxTree, GSProject>) maker.getGoal("taxtree"),
				(ObjectGoal<KMerStoreWrapper, GSProject>) maker.getGoal("updateddb"), fastaTransformGoal);

		//matchGoal.make();
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
