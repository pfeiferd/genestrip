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
import java.nio.charset.StandardCharsets;

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

/**
 * Reads an NCBI assembly summary file (from Genbank or RefSeq) and extracts the assembly entries
 * relevant to a set of taxonomy nodes, filtering by assembly quality and reference status.
 */
public class AssemblySummaryReader {
	/** File name of the RefSeq assembly summary file. */
	public static final String ASSEMLY_SUM_REFSEQ = "assembly_summary_refseq.txt";
	/** File name of the Genbank assembly summary file. */
	public static final String ASSEMLY_SUM_GENBANK = "assembly_summary_genbank.txt";

	/** RefSeq category value marking a reference genome. */
	public static final String REFERENCE_GENOME_CAT = "reference genome";

	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#')
			.setDelimiter('\t').setRecordSeparator('\n').build();

	private final File baseDir;
	private final TaxTree taxTree;
	private final String assFileName;

	/**
	 * Creates a reader for the standard summary file of either Genbank or RefSeq.
	 *
	 * @param baseDir the directory containing the summary file
	 * @param genBank if {@code true}, reads the Genbank summary file, otherwise the RefSeq one
	 * @param taxTree the taxonomy tree used to resolve tax ids
	 * @throws java.io.IOException if the summary file cannot be accessed
	 */
	public AssemblySummaryReader(File baseDir, boolean genBank, TaxTree taxTree) throws IOException {
		this(baseDir, genBank ? ASSEMLY_SUM_GENBANK : ASSEMLY_SUM_REFSEQ, taxTree);
	}

	/**
	 * Creates a reader for the given summary file name.
	 *
	 * @param baseDir the directory containing the summary file
	 * @param assFileName the summary file name
	 * @param taxTree the taxonomy tree used to resolve tax ids
	 * @throws java.io.IOException if the summary file cannot be accessed
	 */
	public AssemblySummaryReader(File baseDir, String assFileName, TaxTree taxTree) throws IOException {
		this.baseDir = baseDir;
		this.taxTree = taxTree;
		this.assFileName = assFileName;
	}

	/**
	 * Parses the summary file and groups the assembly entries by taxonomy node.
	 *
	 * @param filter                      only entries whose node is contained here are kept; {@code null} keeps all
	 * @param fastaQualities              only entries with a quality in this list are kept; {@code null} keeps all
	 * @param referenceOnly               if {@code true}, keeps only reference genomes
	 * @param useSpeciesTaxidForFiltering if {@code true}, resolves and filters by the species tax id
	 *                                    rather than the entry's own tax id
	 * @param totalEntries                if non-{@code null}, its first element receives the total number
	 *                                    of processed records
	 * @return a map from taxonomy node to its relevant assembly entries
	 * @throws java.io.IOException if the summary file cannot be read
	 */
	public Map<TaxIdNode, List<AssemblyEntry>> getRelevantEntries(final Set<TaxIdNode> filter,
																  final List<AssemblyQuality> fastaQualities, final boolean referenceOnly, final boolean useSpeciesTaxidForFiltering, final int[] totalEntries) throws IOException {
		Map<TaxIdNode, List<AssemblyEntry>> result = new HashMap<TaxIdNode, List<AssemblyEntry>>();
		int counter = 0;
		try (CSVParser parser = FORMAT
				.parse(new InputStreamReader(StreamProvider.getInputStreamForFile(new File(baseDir, assFileName)), StandardCharsets.UTF_8))) {
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

	/**
	 * A single assembly entry from the summary file, carrying its tax ids, download location, quality
	 * and derived genomic fasta file name.
	 */
	public static class AssemblyEntry {
		private final String taxid;
		private final String ftpURL;
		private final String fileName;
		private final AssemblyQuality quality;
		private final URL url;
		private final boolean isReference;
		private	final String speciesTaxid;

		/**
		 * Creates an entry and derives the genomic fasta file name from either the ftp directory URL
		 * (appending {@code _genomic.fna.gz}) or the last path segment of the given {@code url}.
		 *
		 * @param taxid the entry's tax id
		 * @param ftpURL the ftp directory URL of the assembly (may be {@code null})
		 * @param quality the assembly quality level
		 * @param url the download URL used when {@code ftpURL} is {@code null}
		 * @param isReference whether this entry is a reference genome
		 * @param speciesTaxid the species-level tax id
		 */
		public AssemblyEntry(String taxid, String ftpURL, AssemblyQuality quality, URL url, boolean isReference, String speciesTaxid) {
			this.taxid = taxid;
			this.ftpURL = ftpURL;
			this.quality = quality;
			if (ftpURL != null) {
				String baseName;
				int pos = ftpURL.lastIndexOf('/');
				if (pos == ftpURL.length() - 1) {
					pos = ftpURL.lastIndexOf('/', pos - 1);
					baseName = ftpURL.substring(pos + 1, ftpURL.length() - 1);
				}
				else {
					baseName = ftpURL.substring(pos + 1);
				}
				this.fileName = baseName + "_genomic.fna.gz";
			} else {
				String path = url.getFile();
				this.fileName = path.substring(path.lastIndexOf('/') + 1);
			}
			this.url = url;
			this.isReference = isReference;
			this.speciesTaxid = speciesTaxid;
		}

		/**
		 * Returns the species-level tax id of this assembly entry.
		 *
		 * @return the species-level tax id
		 */
		public String getSpeciesTaxid() {
			return speciesTaxid;
		}

		/**
		 * Returns the tax id of this assembly entry.
		 *
		 * @return the entry's tax id
		 */
		public String getTaxid() {
			return taxid;
		}

		/**
		 * Indicates whether this entry is a reference genome.
		 *
		 * @return whether this entry is a reference genome
		 */
		public boolean isReference() {
			return isReference;
		}

		/**
		 * Returns the FTP directory URL of the assembly.
		 *
		 * @return the ftp directory URL of the assembly (may be {@code null})
		 */
		public String getFtpURL() {
			return ftpURL;
		}

		/**
		 * Returns the download URL of the assembly.
		 *
		 * @return the download URL of the assembly (may be {@code null})
		 */
		public URL getURL() {
			return url;
		}

		/**
		 * Returns the assembly quality level of this entry.
		 *
		 * @return the assembly quality level
		 */
		public AssemblyQuality getQuality() {
			return quality;
		}

		/**
		 * Returns the derived genomic fasta file name for this assembly entry.
		 *
		 * @return the derived genomic fasta file name
		 */
		public String getFileName() {
			return fileName;
		}
	}

	/**
	 * The assembly quality levels, combining the assembly level (complete, chromosome, scaffold,
	 * contig) with whether the version is the latest.
	 */
	public enum AssemblyQuality {
		/** Additional assembly, treated as the highest quality level. */
		ADDITIONAL,
		/** Complete genome that is the latest version. */
		COMPLETE_LATEST,
		/** Complete genome that is not the latest version. */
		COMPLETE,
		/** Chromosome-level assembly that is the latest version. */
		CHROMOSOME_LATEST,
		/** Chromosome-level assembly that is not the latest version. */
		CHROMOSOME,
		/** Scaffold-level assembly that is the latest version. */
		SCAFFOLD_LATEST,
		/** Scaffold-level assembly that is not the latest version. */
		SCAFFOLD,
		/** Contig-level assembly that is the latest version. */
		CONTIG_LATEST,
		/** Contig-level assembly that is not the latest version. */
		CONTIG,
		/** Latest version with an otherwise unspecified assembly level. */
		LATEST,
		/** No quality level, treated as the lowest quality level. */
		NONE;

		/**
		 * Derives the quality level from the summary file's {@code assembly_level} and
		 * {@code version_status} columns.
		 *
		 * @param complete the {@code assembly_level} column value
		 * @param latest the {@code version_status} column value
		 * @return the derived assembly quality level
		 */
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
