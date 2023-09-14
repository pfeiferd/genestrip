package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;

public class AccessionCollectionGoal extends ObjectGoal<Map<String, TaxIdNode>, GSProject> {
	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setDelimiter('\t').setRecordSeparator('\n')
			.build();

	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final Collection<RefSeqCategory> categories;
	private final File catalogFile;

	@SafeVarargs
	public AccessionCollectionGoal(GSProject project, String name, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			Collection<RefSeqCategory> categories, RefSeqCatalogDownloadGoal catalogGoal, Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, catalogGoal));
		this.taxTreeGoal = taxTreeGoal;
		this.categories = categories;
		catalogFile = catalogGoal.getCatalogFile();
	}

	@Override
	public void makeThis() {
		try {
			Map<String, TaxIdNode> map = new HashMap<String, TaxTree.TaxIdNode>();
			TaxTree taxTree = taxTreeGoal.get();
			Reader in = new InputStreamReader(StreamProvider.getInputStreamForFile(catalogFile));
			Iterable<CSVRecord> records = FORMAT.parse(in);
			for (CSVRecord record : records) {
				String catStr = record.get(3);
				if (containsCategory(catStr)) {
					String taxid = record.get(0);
					if (taxid != null) {
						TaxIdNode node = taxTree.getNodeByTaxId(taxid);
						if (node != null) {
							String accession = record.get(2);
							map.put(accession, node);
						}
					}
				}
			}
			in.close();
			set(map);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

//	protected boolean isIncluded(String taxid, TaxTree taxTree, Set<TaxIdNode> taxNodes) {
//		TaxIdNode node = taxTree.getNodeByTaxId(taxid);
//		return node != null && (taxNodes.isEmpty() || taxNodes.contains(node));
//	}

	protected boolean containsCategory(String catStr) {
		for (RefSeqCategory cat : categories) {
			if (catStr.contains(cat.getDirectory())) {
				return true;
			}
		}
		return false;
	}
}
