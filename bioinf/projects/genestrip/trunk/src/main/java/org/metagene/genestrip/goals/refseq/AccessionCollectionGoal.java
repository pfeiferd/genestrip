package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;

public class AccessionCollectionGoal extends ObjectGoal<List<String>, GSProject> {
	private static final CSVFormat format = CSVFormat.DEFAULT.builder().setDelimiter('\t').setRecordSeparator('\n')
			.build();

	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final Collection<RefSeqCategory> categories;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final File catalogFile;

	@SafeVarargs
	public AccessionCollectionGoal(GSProject project, String name, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			Collection<RefSeqCategory> categories, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			RefSeqCatalogDownloadGoal catalogGoal, Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, taxNodesGoal, catalogGoal));
		this.taxTreeGoal = taxTreeGoal;
		this.taxNodesGoal = taxNodesGoal;
		this.categories = categories;
		catalogFile = catalogGoal.getCatalogFile();
	}

	@Override
	public void makeThis() {
		try {
			List<String> list = new ArrayList<String>();
			TaxTree taxTree = taxTreeGoal.get();
			Set<TaxIdNode> taxNodes = taxNodesGoal.get();
			Reader in = new InputStreamReader(StreamProvider.getInputStreamForFile(catalogFile));
			Iterable<CSVRecord> records = format.parse(in);
			for (CSVRecord record : records) {
				String catStr = record.get(3);
				if (containsCategory(catStr)) {
					String taxid = record.get(0);
					if (isIncluded(taxid, taxTree, taxNodes)) {
						String accession = record.get(2);
						list.add(accession);
					}
				}
			}
			in.close();
			Collections.sort(list);
			set(list);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected boolean isIncluded(String taxid, TaxTree taxTree, Set<TaxIdNode> taxNodes) {
		TaxIdNode node = taxTree.getNodeByTaxId(taxid);
		return node != null && (taxNodes.isEmpty() || taxNodes.contains(node));
	}

	protected boolean containsCategory(String catStr) {
		for (RefSeqCategory cat : categories) {
			if (catStr.contains(cat.name())) {
				return true;
			}
		}
		return false;
	}
}
