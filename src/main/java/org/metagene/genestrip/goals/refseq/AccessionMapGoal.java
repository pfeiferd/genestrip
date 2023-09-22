package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;

public class AccessionMapGoal extends ObjectGoal<AccessionMap, GSProject> {
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final RefSeqCatalogDownloadGoal catalogGoal;
	private final ObjectGoal<Integer, GSProject> accessionMapSizeGoal;
	private final ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal;

	@SafeVarargs
	public AccessionMapGoal(GSProject project, String name, ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal,
			ObjectGoal<TaxTree, GSProject> taxTreeGoal, RefSeqCatalogDownloadGoal catalogGoal,
			RefSeqFnaFilesDownloadGoal downloadGoal, ObjectGoal<Integer, GSProject> accessionMapSizeGoal,
			Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, categoriesGoal, catalogGoal, downloadGoal, accessionMapSizeGoal));
		this.categoriesGoal = categoriesGoal;
		this.taxTreeGoal = taxTreeGoal;
		this.catalogGoal = catalogGoal;
		this.accessionMapSizeGoal = accessionMapSizeGoal;
	}

	@Override
	public void makeThis() {
		AccessionFileProcessor processor = new AccessionFileProcessor(categoriesGoal.get()[1],
				getProject().isUseCompletGenomesOnly()) {
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
		processor.processCatalog(catalogGoal.getCatalogFile());
	}
}
