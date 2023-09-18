package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.util.Collection;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;

public class AccessionMapGoal extends ObjectGoal<AccessionMap, GSProject> {
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final ObjectGoal<Integer, GSProject> accessionMapSizeGoal;
	private final Collection<RefSeqCategory> categories;
	private final File catalogFile;

	@SafeVarargs
	public AccessionMapGoal(GSProject project, String name, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			RefSeqCatalogDownloadGoal catalogGoal, RefSeqFnaFilesDownloadGoal downloadGoal,
			ObjectGoal<Integer, GSProject> accessionMapSizeGoal, Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, catalogGoal, downloadGoal, accessionMapSizeGoal));
		this.taxTreeGoal = taxTreeGoal;
		this.accessionMapSizeGoal = accessionMapSizeGoal;
		this.categories = downloadGoal.getCategories();
		catalogFile = catalogGoal.getCatalogFile();
	}

	@Override
	public void makeThis() {
		AccessionFileProcessor processor = new AccessionFileProcessor(categories) {
			private AccessionMap map = new AccesionMapImpl(accessionMapSizeGoal.get());
			private TaxTree taxTree = taxTreeGoal.get();
			
			@Override
			public void processCatalog(File catalogFile) {
				super.processCatalog(catalogFile);
				map.optimize();
				set(map);
			}
			
			@Override
			protected void handleEntry(byte[] target, int taxIdEnd, int accessionStart, int accessionEnd) {
				TaxIdNode node = taxTree.getNodeByTaxId(target, 0, taxIdEnd);
				if (node != null) {
					map.put(target, accessionStart, accessionEnd, node);
				}
			}
		};
		processor.processCatalog(catalogFile);
	}
}
