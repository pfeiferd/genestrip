package org.metagene.genestrip.match;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.match.FastqKMerMatcher.Result;
import org.metagene.genestrip.match.FastqKMerMatcher.StatsPerTaxid;
import org.metagene.genestrip.tax.TaxTree.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class ResultReporter {
	private final List<TaxIdNode> taxids;

	public ResultReporter(List<TaxIdNode> taxids) {
		this.taxids = taxids;
	}

	public void printStoreInfo(Map<String, Long> counts, PrintStream out) {
		out.println("name;rank;taxid;stored kmers;");
		out.print("TOTAL;");
		out.print(Rank.SUPERKINGDOM);
		out.print(';');
		out.print("1;");
		out.print(counts.get(null));
		out.println(';');
		for (TaxIdNode taxNode : taxids) {
			Long count = counts.get(taxNode.getTaxId());
			if (count == null) {
				count = 0L;
			}
			out.print(taxNode.getName());
			out.print(';');
			out.print(taxNode.getRank());
			out.print(';');
			out.print(taxNode.getTaxId());
			out.print(';');
			out.print(count);
			out.println(';');
		}
	}

	public void print(Result res, PrintStream out) {
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
