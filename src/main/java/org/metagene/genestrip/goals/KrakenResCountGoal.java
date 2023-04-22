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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.kraken.KrakenExecutor;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CountingDigitTrie;
import org.metagene.genestrip.util.StreamProvider;

public class KrakenResCountGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;

	@SafeVarargs
	public KrakenResCountGoal(GSProject project, String name, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			Goal<GSProject>... deps) {
		super(project, name, project.getTaxCountsFile(name), ArraysUtil.append(deps, taxNodesGoal));
		this.taxNodesGoal = taxNodesGoal;
	}

	@Override
	protected void makeFile(File file) throws IOException {
		CountingDigitTrie kmerCountTrie = new CountingDigitTrie();
		CountingDigitTrie classCountTrie = new CountingDigitTrie();

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
							String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output) {
						if (taxIds == null || taxIds.contains(kmerTaxid)) {
							kmerCountTrie.add(kmerTaxid, hitLength);
							System.out.println(ByteArrayUtil.toString(output));
						}
						if (lineCount != lastLine) {
							lastLine = lineCount;
							if (taxIds == null || taxIds.contains(krakenTaxid)) {
//								System.out.println("Classification: " + ByteArrayUtil.toString(output));
//								System.out.println(krakenTaxid + " " + bps + " " + kmerTaxid);
								classCountTrie.inc(krakenTaxid);
							}
						}
					}
				});
			}
		};
		if (getLogger().isInfoEnabled()) {
			String execLine = krakenExecutor.genExecLine(getProject().getKrakenDB(), getProject().getFastqFile());
			getLogger().info("Run kraken with " + execLine);
		}
		try {
			krakenExecutor.execute2(getProject().getKrakenDB(), getProject().getFastqFile(), null, System.err);

			PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));

			Map<String, Long> map = new HashMap<String, Long>();
			out.println("taxid;kmers");
			kmerCountTrie.collect(map);
			CountingDigitTrie.print(map, out);

			map.clear();
			out.println("taxid;classifications");
			classCountTrie.collect(map);
			CountingDigitTrie.print(map, out);

			out.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Finished kraken");
		}
	}
}
