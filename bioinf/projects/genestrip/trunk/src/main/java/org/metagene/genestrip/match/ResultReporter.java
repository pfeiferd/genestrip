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
