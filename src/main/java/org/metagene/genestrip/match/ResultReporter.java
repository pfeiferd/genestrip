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
package org.metagene.genestrip.match;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class ResultReporter {
	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setDelimiter(';')
			.setRecordSeparator('\n').build();

	private static final DecimalFormat DF = new DecimalFormat("0.00000000", new DecimalFormatSymbols(Locale.US));

	public ResultReporter() {
	}

	public void printStoreInfo(Database database, PrintStream out) {
		Object2LongMap<String> stats = database.getStats();

		out.println("name;rank;taxid;stored kmers;");

		out.print("TOTAL;");
		out.print(Rank.NO_RANK);
		out.print(';');
		out.print("1;");
		out.print(stats.getLong(null));
		out.println(';');

		List<String> sortedTaxIds = new ArrayList<String>(stats.keySet());
		SmallTaxTree taxTree = database.getTaxTree();
		taxTree.sortTaxidsViaTree(sortedTaxIds);

		for (String taxId : sortedTaxIds) {
			if (taxId != null) {
				SmallTaxIdNode taxNode = taxTree.getNodeByTaxId(taxId);
				if (taxNode != null) {
					out.print(taxNode.getName());
					out.print(';');
					out.print(taxNode.getRank());
					out.print(';');
					out.print(taxNode.getTaxId());
					out.print(';');
					out.print(stats.getLong(taxId));
					out.println(';');
				}
			}
		}
	}

	/*
	public static List<CountsPerTaxid> readResultCSV(InputStream in, int k) throws IOException {
		Iterable<CSVRecord> records = FORMAT.parse(new InputStreamReader(in));

		List<CountsPerTaxid> res = new ArrayList<CountsPerTaxid>();
		boolean first = true;
		for (CSVRecord record : records) {
			if (first) {
				first = false;
				continue;
			}
			String taxid = record.get(2);
			String reads = record.get(3);
			String kmers = record.get(4);
			String uniqueKmers = record.get(5);
			String contigs = record.get(6);
			String maxContigLen = record.get(8);

			CountsPerTaxid stats = new CountsPerTaxid(taxid, 0);
			stats.reads = Long.parseLong(reads);
			stats.kmers = Long.parseLong(kmers);
			stats.uniqueKmers = Long.parseLong(uniqueKmers);
			stats.contigs = Integer.parseInt(contigs);
			stats.maxContigLen = Integer.parseInt(maxContigLen) - k + 1;
			res.add(stats);
		}

		return res;
	}
	*/

	public void printMatchResult(MatchingResult res, SmallTaxTree taxTree, PrintStream out) {
		out.print("name;rank;taxid;reads;kmers from reads;kmers;unique kmers;contigs;average contig length;max contig length;normalized kmers;exp. unique kmers;unique kmers / exp.;");
		out.print("db coverage;reads bps;avg read len;reads >= 1 kmer;db kmers;parent tax id;");
		if (res.isWithMaxKMerCounts()) {
			out.print("max kmer counts;");
		}
		out.println();
		out.print("TOTAL;");
		out.print(Rank.SUPERKINGDOM);
		out.print(";1;");
		out.print(res.getTotalReads());
		out.print(";0;");
		out.print(res.getTotalKMers());
		out.print(";0;0;0;0;0;0;0;0;0;0;0;0;0;");
		out.print(res.getDbKMers());
		out.print(";1;");
		if (res.isWithMaxKMerCounts()) {
			short[] counts = res.getTotalMaxCounts();
			if (counts != null) {
				for (int i = 0; i < counts.length; i++) {
					out.print(counts[i]);
					out.print(';');
				}
			}
			out.print(';');
		}
		out.println();

		Map<String, CountsPerTaxid> taxid2Stats = res.getTaxid2Stats();
		List<String> keys = new ArrayList<String>(taxid2Stats.keySet());
		taxTree.sortTaxidsViaTree(keys);

		for (String taxId : keys) {
			CountsPerTaxid stats = taxid2Stats.get(taxId);
			if (stats != null) {
				SmallTaxIdNode taxNode = taxTree.getNodeByTaxId(taxId);
				if (taxNode != null) {
					out.print(taxNode.getName());
					out.print(';');
					out.print(taxNode.getRank());
					out.print(';');
					out.print(taxNode.getTaxId());
					out.print(';');
					out.print(stats.getReads());
					out.print(';');
					out.print(stats.getReadKMers());
					out.print(';');
					out.print(stats.getKMers());
					out.print(';');
					out.print(stats.getUniqueKMers());
					out.print(';');
					out.print(stats.getContigs());
					out.print(';');
					out.print(DF.format(((double) stats.getKMers()) / stats.getContigs() + res.getK() - 1));
					out.print(';');
					out.print(stats.getMaxContigLen() + res.getK() - 1);
					out.print(';');
					out.print(DF.format(stats.getCoverage()));
					out.print(';');
					double exp = stats.getExpectedUniqueKMers();
					out.print(DF.format(exp));
					out.print(';');
					out.print(DF.format(stats.getUniqueKMers() / exp));
					out.print(';');
					out.print(stats.getReadBPs());
					out.print(';');
					out.print(DF.format(stats.getAverageReadLength()));
					out.print(';');
					out.print(stats.getReads1KMer());
					out.print(';');
					out.print(stats.getDbKMers());
					out.print(';');
					out.print(stats.getParentTaxId());
					out.print(';');
					ByteArrayUtil.print(stats.getMaxContigDescriptor(), out);
					out.print(';');
					if (res.isWithMaxKMerCounts()) {
						short[] counts = stats.getMaxKMerCounts();
						if (counts != null) {
							for (int i = 0; i < counts.length; i++) {
								if (i > 0) {
									out.print(';');									
								}
								out.print(counts[i]);
							}
						}
						out.print(';');
					}
					out.println();
				}
			}
		}
	}
}
