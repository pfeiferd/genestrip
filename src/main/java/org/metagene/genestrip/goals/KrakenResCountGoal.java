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
package org.metagene.genestrip.goals;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.kraken.KrakenExecutor;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.DigitTrie;
import org.metagene.genestrip.util.StreamProvider;
import org.metagene.genestrip.util.StringLongDigitTrie;

public class KrakenResCountGoal extends FileListGoal<GSProject> {
	private final File fastq;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;

	@SafeVarargs
	public KrakenResCountGoal(GSProject project, String name, File fastqFile,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, fastqFile, FileType.CSV, false),
				ArraysUtil.append(deps, taxNodesGoal));
		this.taxNodesGoal = taxNodesGoal;
		this.fastq = fastqFile;
	}

	@Override
	protected void makeFile(File file) throws IOException {
		DigitTrie<KrakenResStats> countingTrie = new DigitTrie<KrakenResStats>() {
			protected KrakenResStats createInGet(String taxid) {
				KrakenResStats stats = new KrakenResStats();
				stats.taxid = taxid;
				return stats;
			}
		};

		final Set<String> taxIds;

		if (taxNodesGoal != null) {
			taxIds = new HashSet<String>();
			for (TaxIdNode node : taxNodesGoal.get()) {
				taxIds.add(node.getTaxId());
			}
		} else {
			taxIds = null;
		}

		KrakenExecutor krakenExecutor = new KrakenExecutor(getProject().getConfig().getKrakenBin(),
				getProject().getConfig().getKrakenExecExpr()) {
			@Override
			protected void handleOutputStream(InputStream stream, OutputStream out) throws IOException {
				KrakenResultFastqMerger parser = new KrakenResultFastqMerger(
						getProject().getConfig().getMaxReadSizeBytes());

				parser.process(new BufferedInputStream(stream), null, new KrakenResultFastqMergeListener() {
					private long lastLine = -1;

					@Override
					public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
							String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output,
							StringLongDigitTrie root) {
						if (taxIds == null || taxIds.contains(kmerTaxid)) {
							KrakenResStats stats = countingTrie.get(kmerTaxid, true);
							stats.kmers += hitLength;
						}
						if (lineCount != lastLine) {
							lastLine = lineCount;
							if (taxIds == null || taxIds.contains(krakenTaxid)) {
								KrakenResStats stats = countingTrie.get(krakenTaxid, true);
								stats.reads++;
								if (kmerTaxid != null && kmerTaxid.equals(krakenTaxid)) {
									stats.kmersInMatchingReads += hitLength;
								}
							}
						}
					}
				});
			}
		};
		if (getLogger().isInfoEnabled()) {
			String execLine = krakenExecutor.genExecLine(getProject().getKrakenDB(), fastq, file);
			getLogger().info("Run kraken with " + execLine);
		}
		try {
			if (krakenExecutor.isWithFileForOutput()) {
				throw new IOException("This goal does not work with an outfile as a parameter (like in krakenuniq)");
			}

			krakenExecutor.execute2(getProject().getKrakenDB(), fastq, file, null, System.err);

			PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));

			List<KrakenResStats> list = new ArrayList<KrakenResStats>();
			countingTrie.collect(list);
			print(list, out);

			out.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Finished kraken");
		}
	}

	protected void print(List<KrakenResStats> list, PrintStream out) {
		out.println("taxid;reads;kmers;kmers in matching reads");
		for (KrakenResStats stats : list) {
			out.print(stats.taxid);
			out.print(';');
			out.print(stats.reads);
			out.print(';');
			out.print(stats.kmers);
			out.print(';');
			out.print(stats.kmersInMatchingReads);
			out.println(';');
		}
	}

	protected static class KrakenResStats implements Comparable<KrakenResStats> {
		protected String taxid;

		protected long reads;
		protected long kmers;
		protected long kmersInMatchingReads;

		@Override
		public int compareTo(KrakenResStats o) {
			return taxid.compareTo(o.taxid);
		}
	}
}
