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

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.match.FastqKMerMatcher.Result;
import org.metagene.genestrip.match.FastqKMerMatcher.StatsPerTaxid;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.store.KMerStoreWrapper.StoreStatsPerTaxid;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class ResultReporter {
	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setDelimiter(';')
			.setRecordSeparator('\n').build();

	private final TaxTree taxTree;

	public ResultReporter(TaxTree taxTree) {
		this.taxTree = taxTree;
	}

	public void printStoreInfo(Object2LongMap<String> stats, PrintStream out) {
		out.println("name;rank;taxid;stored kmers;");

		out.print("TOTAL;");
		out.print(Rank.SUPERKINGDOM);
		out.print(';');
		out.print("1;");
		out.print(stats.getLong(null));
		out.println(';');

		List<String> sortedTaxIds = new ArrayList<String>(stats.keySet());
		taxTree.sortTaxidsViaTree(sortedTaxIds);
		
		for (String taxId : sortedTaxIds) {
			if (taxId != null) {
				TaxIdNode taxNode = taxTree.getNodeByTaxId(taxId);
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

	// TODO: Outdated..
	public static List<StoreStatsPerTaxid> readStoreInfoCSV(InputStream in) throws IOException {
		Iterable<CSVRecord> records = FORMAT.parse(new InputStreamReader(in));

		List<StoreStatsPerTaxid> res = new ArrayList<StoreStatsPerTaxid>();
		boolean first = true;
		for (CSVRecord record : records) {
			if (first) {
				first = false;
				continue;
			}
			String taxid = record.get(2);
			String storedKMers = record.get(3);
			String totalKMers = record.get(4);
			String assignedKMers = record.get(6);
			String contigs = record.get(7);
			String maxContigLen = record.get(9);

			StoreStatsPerTaxid stats = new StoreStatsPerTaxid(taxid);
			stats.storedKMers = Long.parseLong(storedKMers);
			stats.totalKMers = Long.parseLong(totalKMers);
			stats.assignedKMers = Long.parseLong(assignedKMers);
			stats.contigs = Long.parseLong(contigs);
			stats.maxContigLen = Long.parseLong(maxContigLen);
			res.add(stats);
		}

		return res;
	}

	public static List<StatsPerTaxid> readResultCSV(InputStream in) throws IOException {
		Iterable<CSVRecord> records = FORMAT.parse(new InputStreamReader(in));

		List<StatsPerTaxid> res = new ArrayList<StatsPerTaxid>();
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

			StatsPerTaxid stats = new StatsPerTaxid(taxid);
			stats.reads = Long.parseLong(reads);
			stats.kmers = Long.parseLong(kmers);
			stats.uniqueKmers = Long.parseLong(uniqueKmers);
			stats.contigs = Integer.parseInt(contigs);
			stats.maxContigLen = Integer.parseInt(maxContigLen);
			res.add(stats);
		}

		return res;
	}

	public void printMatchResult(Result res, PrintStream out, KMerStoreWrapper wrapper) {
		Map<String, StatsPerTaxid> taxid2Stats = res.getTaxid2Stats();
		UniqueKMerEstimator estimator = wrapper == null ? null : new UniqueKMerEstimator(wrapper);
		if (estimator != null) {
			estimator.setTotalKMers(res.getTotalKMers());
		}

		out.print("name;rank;taxid;reads;kmers;unique kmers;contigs;average contig length;max contig length;");
		if (estimator != null) {
			out.print(
					"normalized kmers; exp. unique kmers; unique kmers / exp.; quality prediction;");
		}
		if (res.isWithMaxKMerCounts()) {
			out.print("max kmer counts;");
		}
		out.println();
		out.print("TOTAL;");
		out.print(Rank.SUPERKINGDOM);
		out.print(';');
		out.print("1;");
		out.print(res.getTotalReads());
		out.print(';');
		out.print(res.getTotalKMers());
		out.print(";0;0;0;0;");
		if (estimator != null) {
			out.print("0;0;0;0;");
		}
		if (res.isWithMaxKMerCounts()) {
			short[] counts = res.getTotalMaxCounts();
			if (counts != null) {
				for (int i = 0; i < counts.length; i++) {
					out.print(counts[i]);
					out.print('|');
				}
			}
			out.print(';');
		}
		out.println();

		List<String> sortedTaxIds = new ArrayList<String>(taxid2Stats.keySet());
		taxTree.sortTaxidsViaTree(sortedTaxIds);
		for (String taxId : sortedTaxIds) {
			StatsPerTaxid stats = taxid2Stats.get(taxId);
			if (stats != null) {
				TaxIdNode taxNode = taxTree.getNodeByTaxId(taxId);
				if (taxNode != null) {
					out.print(taxNode.getName());
					out.print(';');
					out.print(taxNode.getRank());
					out.print(';');
					out.print(taxNode.getTaxId());
					out.print(';');
					out.print(stats.getReads());
					out.print(';');
					out.print(stats.getKMers());
					out.print(';');
					out.print(stats.getUniqueKMers());
					out.print(';');
					out.print(stats.getContigs());
					out.print(';');
					out.print(((double) stats.getKMers()) / stats.getContigs());
					out.print(';');
					out.print(stats.getMaxContigLen());
					out.print(';');
					if (estimator != null) {
						double normalizedKMers = estimator.getNormalizedKMers(stats);
						out.print(normalizedKMers);
						out.print(';');
						double expUnique = estimator.getExpectedUniqueKMers(stats);
						out.print(expUnique);
						out.print(';');
						double uniqueExpRatio = stats.getUniqueKMers() / expUnique;
						out.print(uniqueExpRatio);
						out.print(';');
						out.print(normalizedKMers * uniqueExpRatio);
						out.print(';');
					}
					if (res.isWithMaxKMerCounts()) {
						short[] counts = stats.getMaxKMerCounts();
						if (counts != null) {
							for (int i = 0; i < counts.length; i++) {
								out.print(counts[i]);
								out.print('|');
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
