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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.CountingDigitTrie;
import org.metagene.genestrip.util.StreamProvider;

public class KrakenFastqFileGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final FileGoal<GSProject> krakenOutGoal;
	private final KMerFastqGoal kmerFastqGoal;
	
	private long lastLineCount;
	private long readCount;

	@SafeVarargs
	public KrakenFastqFileGoal(GSProject project, String name, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			FileGoal<GSProject> krakenOutGoal, KMerFastqGoal kmerFastqGoal, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.FASTQ), ArraysUtil.append(deps, taxNodesGoal, krakenOutGoal, kmerFastqGoal));
		this.taxNodesGoal = taxNodesGoal;
		this.krakenOutGoal = krakenOutGoal;
		this.kmerFastqGoal = kmerFastqGoal;
	}

	@Override
	protected void makeFile(File fastq) {
		try {
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Writing file " + fastq);
			}
			PrintStream printStream = new PrintStream(StreamProvider.getOutputStreamForFile(fastq));

			KrakenResultFastqMergeListener filter = KrakenResultFastqMergeListener
					.createFilterByTaxIdNodes(taxNodesGoal.get(), new KrakenResultFastqMergeListener() {
						@Override
						public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read,
								byte[] readProbs, String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength,
								byte[] output, CountingDigitTrie root) {
							if (lastLineCount == lineCount) {
								throw new IllegalStateException("too many k-mers for read");
							}
							lastLineCount = lineCount;
							readCount++;
							ByteArrayUtil.print(readDescriptor, printStream);
							printStream.print(":taxid:");
							printStream.println(kmerTaxid);
							ByteArrayUtil.print(read, printStream);
							printStream.println();
							printStream.println("+");
							ByteArrayUtil.print(readProbs, printStream);
							printStream.println();
						}
					});
			
			lastLineCount = 0;
			readCount = 0;
			KrakenResultFastqMerger krakenKMerFastqMerger = new KrakenResultFastqMerger(
					getProject().getConfig().getMaxReadSizeBytes());

			for (int i = 0; i < krakenOutGoal.getFiles().size(); i++) {
				if (getLogger().isInfoEnabled()) {
					getLogger().info("Reading file " + krakenOutGoal.getFiles().get(i));
					getLogger().info("Reading file " + kmerFastqGoal.getFiles().get(i));
				}
				
				InputStream stream1 = StreamProvider.getInputStreamForFile(krakenOutGoal.getFiles().get(i));
				InputStream stream2 = StreamProvider.getInputStreamForFile(kmerFastqGoal.getFiles().get(i));
				krakenKMerFastqMerger.process(stream1, stream2, filter);
				stream1.close();
				stream2.close();
			}

			printStream.close();
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Written file " + fastq);
				getLogger().info("Total added reads: " + readCount);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

}
