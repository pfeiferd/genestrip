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
package org.metagene.genestrip.genbank;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class AssemblySummaryReader {
	public static final String ASSEMLY_SUM_REFSEQ = "assembly_summary_refseq.txt";
	public static final String ASSEMLY_SUM_GENBANK = "assembly_summary_genbank.txt";

	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#')
			.setDelimiter('\t').setRecordSeparator('\n').build();

	private final File baseDir;
	private final TaxTree taxTree;
	private final String assFileName;

	public AssemblySummaryReader(File baseDir, boolean genBank, TaxTree taxTree) throws IOException {
		this(baseDir, genBank ? ASSEMLY_SUM_GENBANK : ASSEMLY_SUM_REFSEQ, taxTree);
	}

	public AssemblySummaryReader(File baseDir, String assFileName, TaxTree taxTree) throws IOException {
		this.baseDir = baseDir;
		this.taxTree = taxTree;
		this.assFileName = assFileName;
	}

	public Map<TaxIdNode, List<FTPEntryWithQuality>> getRelevantEntries(Set<TaxIdNode> filter,
			List<FTPEntryQuality> minQualities, int[] totalEntries) throws IOException {
		Map<TaxIdNode, List<FTPEntryWithQuality>> result = new HashMap<TaxIdNode, List<FTPEntryWithQuality>>();
		int counter = 0;

		try (CSVParser parser = FORMAT
				.parse(new InputStreamReader(StreamProvider.getInputStreamForFile(new File(baseDir, assFileName))))) {
			for (CSVRecord record : parser) {
				// Prevent inconsistent entries (they gotta have at least 20 columns)
				if (record.size() >= 20) {
					String taxid = record.get(5);
					String latest = record.get(10);
					String complete = record.get(11);
					String ftp = record.get(19);

					TaxIdNode node = taxTree.getNodeByTaxId(taxid);
					if (node != null && (filter == null || filter.contains(node))) {
						FTPEntryQuality quality = FTPEntryQuality.fromString(complete, latest);
						if (!FTPEntryQuality.NONE.equals(quality)) {
							List<FTPEntryWithQuality> entry = result.get(node);
							if (entry == null) {
								entry = new ArrayList<FTPEntryWithQuality>();
								result.put(node, entry);
							}
							FTPEntryWithQuality ewq = new FTPEntryWithQuality(ftp, quality, null);
							entry.add(ewq);
						}
					}
					counter++;
				}
			}
		}

		if (minQualities != null) {
			for (List<FTPEntryWithQuality> values : result.values()) {
				boolean found = false;
				// Try the qualities in given order:
				for (FTPEntryQuality q : minQualities) {
					for (FTPEntryWithQuality entry : values) {
						if (!entry.getQuality().below(q)) {
							// At least one entry matches the required quality q.
							found = true;
							break;
						}
					}
					if (found) {
						// Remove all entries that don't match the quality q
						Iterator<FTPEntryWithQuality> it = values.iterator();
						while (it.hasNext()) {
							if (it.next().getQuality().below(q)) {
								it.remove();
							}
						}

						break;
					}
				}
				// If not matches were found, then the required qualities are exhausted and we
				// got to delete all entries.
				if (!found) {
					values.clear();
				}
			}
		}
		if (totalEntries != null) {
			totalEntries[0] = counter;
		}
		return result;
	}

	public static class FTPEntryWithQuality {
		private final String ftpURL;
		private final String fileName;
		private final FTPEntryQuality quality;
		private final URL url;

		public FTPEntryWithQuality(String ftpURL, FTPEntryQuality quality, URL url) {
			this.ftpURL = ftpURL;
			this.quality = quality;
			if (ftpURL != null) {
				this.fileName = ftpURL.substring(ftpURL.lastIndexOf('/') + 1) + "_genomic.fna.gz";
			} else {
				String path = url.getFile();
				this.fileName = path.substring(path.lastIndexOf('/') + 1);
			}
			this.url = url;
		}

		public String getFtpURL() {
			return ftpURL;
		}

		public URL getURL() {
			return url;
		}

		public FTPEntryQuality getQuality() {
			return quality;
		}

		public String getFileName() {
			return fileName;
		}
	}

	public enum FTPEntryQuality {
		ADDITIONAL, COMPLETE_LATEST, COMPLETE, CHROMOSOME_LATEST, CHROMOSOME, CONTIG_LATEST, CONTIG, LATEST, NONE;

		public boolean below(FTPEntryQuality q) {
			return this.ordinal() > q.ordinal();
		}

		public static FTPEntryQuality fromString(String complete, String latest) {
			if ("Complete Genome".equals(complete)) {
				if ("latest".equals(latest)) {
					return COMPLETE_LATEST;
				} else {
					return COMPLETE;
				}
			} else if ("Chromosome".equals(complete)) {
				if ("latest".equals(latest)) {
					return CHROMOSOME_LATEST;
				} else {
					return CHROMOSOME;
				}
			} else if ("Contig".equals(complete)) {
				if ("latest".equals(latest)) {
					return CONTIG_LATEST;
				} else {
					return CONTIG;
				}
			} else {
				if ("latest".equals(latest)) {
					return LATEST;
				} else {
					return NONE;
				}
			}
		}
	}
}
