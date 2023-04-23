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
import java.util.HashSet;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.KMerTrie;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;

public class KMerTrieFileGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final KrakenOutGoal krakenOutGoal;
	private final KMerFastqGoal kmerFastqGoal;

	@SafeVarargs
	public KMerTrieFileGoal(GSProject project, String name, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			KrakenOutGoal krakenOutGoal, KMerFastqGoal kmerFastqGoal, Goal<GSProject>... deps) {
		super(project, name, project.getTrieFile(), ArraysUtil.append(deps, taxNodesGoal, kmerFastqGoal, krakenOutGoal));
		this.taxNodesGoal = taxNodesGoal;
		this.krakenOutGoal = krakenOutGoal;
		this.kmerFastqGoal = kmerFastqGoal;
	}
	
	@Override
	protected void makeFile(File trieFile) {
		try {
			KMerTrie<String> trie = new KMerTrie<String>(1, getProject().getkMserSize(), false);

			Set<String> taxIds = new HashSet<String>();
			for (TaxIdNode node : taxNodesGoal.get()) {
				taxIds.add(node.getTaxId());
			}

			KrakenResultFastqMergeListener filter = KrakenResultFastqMergeListener.createFilterByTaxIds(taxIds,
					fillKMerTrie(trie, null));
			KrakenResultFastqMerger krakenKMerFastqMerger = new KrakenResultFastqMerger(
					getProject().getConfig().getMaxReadSizeBytes());

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Reading file " + krakenOutGoal.getOutputFile());
				getLogger().info("Reading file " + kmerFastqGoal.getOutputFile());
			}
			
			InputStream stream1 = StreamProvider.getInputStreamForFile(krakenOutGoal.getOutputFile());
			InputStream stream2 = StreamProvider.getInputStreamForFile(kmerFastqGoal.getOutputFile());			
			krakenKMerFastqMerger.process(stream1, stream2, filter);
			stream1.close();
			stream2.close();

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Trie entries: " + trie.getEntries());
				getLogger().info("Saving File " + trieFile);
			}
			trie.compress();
			trie.save(trieFile);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File saved " + trieFile);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	
	private KrakenResultFastqMergeListener fillKMerTrie(KMerTrie<String> trie,
			KrakenResultFastqMergeListener delegate) {
		return new KrakenResultFastqMergeListener() {
			@Override
			public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
					String krakenTaxid, int bps, int pos, String kmerTaxid, int hitLength, byte[] output) {
				trie.put(read, 0, kmerTaxid, false);
				if (lineCount % 10000 == 0) {
					if (getLogger().isInfoEnabled()) {
						getLogger().info("Trie entries:" + trie.getEntries());
						getLogger().info("Trie put ratio:" + ((double) trie.getEntries() / lineCount));
					}
				}
				if (delegate != null) {
					delegate.newTaxIdForRead(lineCount, readDescriptor, read, readProbs, krakenTaxid, bps, pos,
							kmerTaxid, hitLength, output);
				}
			}
		};
	}
}
