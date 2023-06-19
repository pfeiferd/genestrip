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
package org.metagene.genestrip.tax;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryQuality;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryWithQuality;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.StreamProvider;

public class AdditionalFastaInfoReader {
	public static final String ADDITIONAL_INFO_FILE = "additional.txt";

	private final File baseDir;
	private final TaxTree taxTree;
	private final String assFileName;

	public AdditionalFastaInfoReader(File baseDir, TaxTree taxTree) throws IOException {
		this(baseDir, ADDITIONAL_INFO_FILE, taxTree);
	}

	public AdditionalFastaInfoReader(File baseDir, String assFileName, TaxTree taxTree) throws IOException {
		this.baseDir = baseDir;
		this.taxTree = taxTree;
		this.assFileName = assFileName;
	}

	public void addRelevantEntries(Map<TaxIdNode, List<FTPEntryWithQuality>> result, Set<TaxIdNode> filter)
			throws IOException {
		addRelevantEntries(result, filter, null);
	}

	public void addRelevantEntries(Map<TaxIdNode, List<FTPEntryWithQuality>> result, Set<TaxIdNode> filter,
			int[] totalEntries) throws IOException {
		File additonalEntryFile = new File(baseDir, assFileName);

		if (additonalEntryFile.exists()) {
			Reader in = new InputStreamReader(StreamProvider.getInputStreamForFile(additonalEntryFile));
			CSVFormat format = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#').setDelimiter('\t')
					.setRecordSeparator('\n').build();
			Iterable<CSVRecord> records = format.parse(in);

			int counter = 0;
			for (CSVRecord record : records) {
				String taxid = record.get(0);
				String url = record.get(1);
				String latest = null;
				if (record.size() > 2) {
					latest = record.get(2);
				}
				String complete = null;
				if (record.size() > 3) {
					complete = record.get(3);
				}

				TaxIdNode node = taxTree.getNodeByTaxId(taxid);
				if (node != null && (filter == null || filter.contains(node))) {
					FTPEntryQuality quality = FTPEntryQuality.fromString(complete, latest);
					if (quality == null || FTPEntryQuality.NONE.equals(quality)) {
						quality = FTPEntryQuality.ADDITIONAL;
					}
					if (!FTPEntryQuality.NONE.equals(quality)) {
						List<FTPEntryWithQuality> entry = result.get(node);
						if (entry == null) {
							entry = new ArrayList<FTPEntryWithQuality>();
							result.put(node, entry);
						}
						FTPEntryWithQuality ewq = new FTPEntryWithQuality(null, quality, new URL(url));
						entry.add(ewq);
					}
				}
				counter++;
			}
			if (totalEntries != null) {
				totalEntries[0] = counter;
			}
		}
	}
}
