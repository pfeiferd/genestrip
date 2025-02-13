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

	public static final String REFERENCE_GENOME_CAT = "reference genome";

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

	public Map<TaxIdNode, List<AssemblyEntry>> getRelevantEntries(final Set<TaxIdNode> filter,
																  final List<AssemblyQuality> fastaQualities, final boolean referenceOnly, final boolean useSpeciesTaxidForFiltering, final int[] totalEntries) throws IOException {
		Map<TaxIdNode, List<AssemblyEntry>> result = new HashMap<TaxIdNode, List<AssemblyEntry>>();
		int counter = 0;
		try (CSVParser parser = FORMAT
				.parse(new InputStreamReader(StreamProvider.getInputStreamForFile(new File(baseDir, assFileName))))) {
			for (CSVRecord record : parser) {
				// Prevent inconsistent entries (they gotta have at least 20 columns)
				if (record.size() >= 20) {
					String refgen = record.get(4);
					String taxid = record.get(5);
					String speciesTaxid = record.get(6);
					String latest = record.get(10);
					String complete = record.get(11);
					String ftp = record.get(19);

					TaxIdNode node = taxTree.getNodeByTaxId(useSpeciesTaxidForFiltering ? speciesTaxid : taxid);
					if (node != null && (filter == null || filter.contains(node))) {
						AssemblyQuality quality = AssemblyQuality.fromString(complete, latest);
						if (fastaQualities == null || fastaQualities.contains(quality)) {
							boolean isRefGen = REFERENCE_GENOME_CAT.equals(refgen);
							if (!referenceOnly || isRefGen) {
								List<AssemblyEntry> entry = result.get(node);
								if (entry == null) {
									entry = new ArrayList<AssemblyEntry>();
									result.put(node, entry);
								}
								AssemblyEntry ewq = new AssemblyEntry(taxid, ftp, quality, null, isRefGen, speciesTaxid);
								entry.add(ewq);
							}
						}
					}
					counter++;
				}
			}
		}
		if (totalEntries != null) {
			totalEntries[0] = counter;
		}
		return result;
	}

	public static class AssemblyEntry {
		private final String taxid;
		private final String ftpURL;
		private final String fileName;
		private final AssemblyQuality quality;
		private final URL url;
		private final boolean isReference;
		private	final String speciesTaxid;

		public AssemblyEntry(String taxid, String ftpURL, AssemblyQuality quality, URL url, boolean isReference, String speciesTaxid) {
			this.taxid = taxid;
			this.ftpURL = ftpURL;
			this.quality = quality;
			if (ftpURL != null) {
				this.fileName = ftpURL.substring(ftpURL.lastIndexOf('/') + 1) + "_genomic.fna.gz";
			} else {
				String path = url.getFile();
				this.fileName = path.substring(path.lastIndexOf('/') + 1);
			}
			this.url = url;
			this.isReference = isReference;
			this.speciesTaxid = speciesTaxid;
		}

		public String getSpeciesTaxid() {
			return speciesTaxid;
		}

		public String getTaxid() {
			return taxid;
		}

		public boolean isReference() {
			return isReference;
		}

		public String getFtpURL() {
			return ftpURL;
		}

		public URL getURL() {
			return url;
		}

		public AssemblyQuality getQuality() {
			return quality;
		}

		public String getFileName() {
			return fileName;
		}
	}

	public enum AssemblyQuality {
		ADDITIONAL, COMPLETE_LATEST, COMPLETE, CHROMOSOME_LATEST, CHROMOSOME, SCAFFOLD_LATEST, SCAFFOLD, CONTIG_LATEST, CONTIG, LATEST, NONE;

		public static AssemblyQuality fromString(String complete, String latest) {
			boolean l = "latest".equals(latest);
			switch (complete) {
				case "Complete Genome":
					return l ? COMPLETE_LATEST: COMPLETE;
				case "Chromosome":
					return l ? CHROMOSOME_LATEST: CHROMOSOME;
				case "Scaffold":
					return l ? SCAFFOLD_LATEST: SCAFFOLD;
				case "Contig":
					return l ? CONTIG_LATEST: CONTIG;
				default:
					return l ? LATEST: NONE;
			}
		}
	}
}
