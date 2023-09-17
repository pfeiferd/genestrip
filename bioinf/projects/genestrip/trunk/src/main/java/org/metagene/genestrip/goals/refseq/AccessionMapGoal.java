package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.BufferedLineReader;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.StreamProvider;
import org.metagene.genestrip.util.StreamProvider.ByteCountingInputStreamAccess;

public class AccessionMapGoal extends ObjectGoal<Map<String, TaxIdNode>, GSProject> {
	protected final Log logger = LogFactory.getLog(getClass());

	private final long recordLogCycle = 1000 * 1000 * 10;

	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final Collection<RefSeqCategory> categories;
	private final File catalogFile;

	@SafeVarargs
	public AccessionMapGoal(GSProject project, String name, ObjectGoal<TaxTree, GSProject> taxTreeGoal, RefSeqCatalogDownloadGoal catalogGoal,
			RefSeqFnaFilesDownloadGoal downloadGoal, Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, catalogGoal, downloadGoal));
		this.taxTreeGoal = taxTreeGoal;
		this.categories = downloadGoal.getCategories();
		catalogFile = catalogGoal.getCatalogFile();
	}

	@Override
	public void makeThis() {
		try {
			ByteCountingInputStreamAccess byteCountAccess = StreamProvider
					.getByteCountingInputStreamForFile(catalogFile, false);

			long totalCatSize = Files.size(catalogFile.toPath());
			byte[] target = new byte[2000];

			TaxTree taxTree = taxTreeGoal.get();

			Map<String, TaxIdNode> map = new HashMap<String, TaxTree.TaxIdNode>();

			// The file is huge, apache csv reader would be too slow and burn too many
			// string.
			// Therefore manual coding for parsing and
			int size;
			BufferedLineReader reader = new BufferedLineReader(byteCountAccess.getInputStream());

			long recordCounter = 0;
			long startTime = System.currentTimeMillis();
			while ((size = reader.nextLine(target)) > 0) {
				int pos1 = ByteArrayUtil.indexOf(target, 0, size, '\t');
				int pos2 = ByteArrayUtil.indexOf(target, pos1 + 1, size, '\t');
				int pos3 = ByteArrayUtil.indexOf(target, pos2 + 1, size, '\t');
				int pos4 = ByteArrayUtil.indexOf(target, pos3 + 1, size, '\t');
				if (containsCategory(target, pos3 + 1, pos4)) {
					TaxIdNode node = taxTree.getNodeByTaxId(target, 0, pos1);
					if (node != null) {
						String accession = new String(target, pos2 + 1, pos3 - pos2 - 1);
						map.put(accession, node);					
					}
				}

				if (recordCounter % recordLogCycle == 0) {
					if (logger.isInfoEnabled()) {
						double ratio = byteCountAccess.getBytesRead() / (double) totalCatSize;
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

			reader.close();
			set(map);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected boolean containsCategory(byte[] outerArray, int start, int end) {
		for (RefSeqCategory cat : categories) {
			if (ByteArrayUtil.indexOf(outerArray, start, end, cat.getDirectory()) != -1) {
				return true;
			}
		}
		return false;
	}
}
