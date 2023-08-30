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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.match.FastqKMerMatcher.Result;
import org.metagene.genestrip.match.FastqKMerMatcher.StatsPerTaxid;
import org.metagene.genestrip.store.KMerStoreWrapper.StoreStatsPerTaxid;
import org.metagene.genestrip.tax.TaxTree.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class ResultReporter {
	private static final CSVFormat format = CSVFormat.DEFAULT.builder().setQuote(null).setDelimiter(';')
			.setRecordSeparator('\n').build();

	private final List<TaxIdNode> taxids;

	public ResultReporter(List<TaxIdNode> taxids) {
		this.taxids = taxids;
	}

	public void printStoreInfo(List<StoreStatsPerTaxid> statsList, PrintStream out) {
		out.println(
				"name;rank;taxid;stored kmers;total kmers;stored ratio;assigned kmers;contigs;average contig length;max contig length;");

		Map<String, StoreStatsPerTaxid> map = new HashMap<String, StoreStatsPerTaxid>();
		for (StoreStatsPerTaxid stats : statsList) {
			map.put(stats.getTaxid(), stats);
		}

		StoreStatsPerTaxid stats = map.get(null);
		out.print("TOTAL;");
		out.print(Rank.SUPERKINGDOM);
		out.print(';');
		out.print("1;");
		out.print(stats.getStoredKMers());
		out.print(';');
		out.print(stats.getTotalKMers());
		out.print(';');
		out.print(((double) stats.getStoredKMers()) / stats.getTotalKMers());
		out.print(';');
		out.print(0);
		out.print(';');
		out.print(0);
		out.print(';');
		out.print('-');
		out.print(';');
		out.print(0);
		out.println(';');

		for (TaxIdNode taxNode : taxids) {
			stats = map.get(taxNode.getTaxId());
			if (stats != null) {
				out.print(taxNode.getName());
				out.print(';');
				out.print(taxNode.getRank());
				out.print(';');
				out.print(taxNode.getTaxId());
				out.print(';');
				out.print(stats.getStoredKMers());
				out.print(';');
				out.print(stats.getTotalKMers());
				out.print(';');
				out.print(((double) stats.getStoredKMers()) / stats.getTotalKMers());
				out.print(';');
				out.print(stats.getAssignedKMers());
				out.print(';');
				out.print(stats.getContigs());
				out.print(';');
				out.print(((double) stats.getStoredKMers()) / stats.getContigs());
				out.print(';');
				out.print(stats.getMaxContigLen());
				out.println(';');
//			} else {
//				out.print(0);
//				out.print(';');
//				out.print(0);
//				out.print(';');
//				out.print(0);
//				out.print(';');
//				out.print('-');
//				out.print(';');
//				out.print('0');
//				out.println(';');
			}
		}
	}

	public static List<StoreStatsPerTaxid> readStoreInfoCSV(InputStream in) throws IOException {
		Iterable<CSVRecord> records = format.parse(new InputStreamReader(in));

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
		Iterable<CSVRecord> records = format.parse(new InputStreamReader(in));

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

	public void printMatchResult(Result res, PrintStream out) {
		Map<String, StatsPerTaxid> taxid2Stats = new HashMap<String, StatsPerTaxid>();
		for (StatsPerTaxid stats : res.getStats()) {
			taxid2Stats.put(stats.getTaxid(), stats);
		}

		out.println("name;rank;taxid;reads;kmers;unique kmers;contigs;average contig length;max contig length;");
		out.print("TOTAL;");
		out.print(Rank.SUPERKINGDOM);
		out.print(';');
		out.print("1;");
		out.print(res.getTotalReads());
		out.print(';');
		out.print(res.getTotalKMers());
		out.println(";0;0;0;0;");

		for (TaxIdNode taxNode : taxids) {
			StatsPerTaxid stats = taxid2Stats.get(taxNode.getTaxId());
			if (stats != null) {
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
				out.println(';');
			}
		}
	}
}
