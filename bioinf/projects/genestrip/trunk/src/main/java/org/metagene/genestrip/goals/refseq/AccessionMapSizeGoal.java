package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.util.Collection;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.util.ArraysUtil;

public class AccessionMapSizeGoal extends ObjectGoal<Integer, GSProject> {
	private final Collection<RefSeqCategory> categories;
	private final RefSeqCatalogDownloadGoal catalogGoal;

	@SafeVarargs
	public AccessionMapSizeGoal(GSProject project, String name, RefSeqCatalogDownloadGoal catalogGoal,
			RefSeqFnaFilesDownloadGoal downloadGoal, Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, catalogGoal, downloadGoal));
		this.catalogGoal = catalogGoal;
		this.categories = downloadGoal.getCategories();
	}

	@Override
	public void makeThis() {
		File catalogFile = catalogGoal.getCatalogFile();
		AccessionFileProcessor processor = new AccessionFileProcessor(categories) {
			private int counter = 0;
			
			@Override
			public void processCatalog(File catalogFile) {
				super.processCatalog(catalogFile);
				set(counter);
				
				if (getLogger().isInfoEnabled()) {
					getLogger().info("Map size determined: " + counter);
				}
			}
			
			@Override
			protected void handleEntry(byte[] target, int taxIdEnd, int accessionStart, int accessionEnd) {
				counter++;
			}
		};
		processor.processCatalog(catalogFile);
	}
}
