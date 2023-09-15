package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.nio.file.Files;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;
import org.metagene.genestrip.util.StreamProvider.ByteCountingInputStreamAccess;

public class AccessionCollectionGoal extends ObjectGoal<Map<String, TaxIdNode>, GSProject> {
	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setDelimiter('\t').setRecordSeparator('\n')
			.build();

	protected final Log logger = LogFactory.getLog(getClass());

	private final long recordLogCycle = 1000 * 1000 * 10;

	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final Collection<RefSeqCategory> categories;
	private final File catalogFile;

	@SafeVarargs
	public AccessionCollectionGoal(GSProject project, String name, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			RefSeqCatalogDownloadGoal catalogGoal, RefSeqFnaFilesDownloadGoal downloadGoal, Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, taxTreeGoal, catalogGoal, downloadGoal));
		this.taxTreeGoal = taxTreeGoal;
		this.categories = downloadGoal.getCategories();
		catalogFile = catalogGoal.getCatalogFile();
	}

	@Override
	public void makeThis() {
		try {
			Map<String, TaxIdNode> map = new HashMap<String, TaxTree.TaxIdNode>();
			TaxTree taxTree = taxTreeGoal.get();
			ByteCountingInputStreamAccess byteCountAccess = StreamProvider
					.getByteCountingInputStreamForFile(catalogFile, false);

			long totalCatSize = Files.size(catalogFile.toPath());

			Reader in = new InputStreamReader(byteCountAccess.getInputStream());
			Iterable<CSVRecord> records = FORMAT.parse(in);
			long recordCounter = 0;
			long startTime = System.currentTimeMillis();
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
				if (recordCounter % recordLogCycle == 0) {
					if (logger.isInfoEnabled()) {
						double ratio =  byteCountAccess.getBytesRead() / (double) totalCatSize;
						long stopTime = System.currentTimeMillis();

						double diff = (stopTime - startTime);
						double totalTime = diff / ratio;
						double totalHours = totalTime / 1000 / 60 / 60;

						logger.info("Elapsed hours:" + diff / 1000 / 60 / 60);
						logger.info("Estimated total hours:" + totalHours);
						logger.info("Records processed: " + recordCounter);
					}
				}
				recordCounter++;
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
