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
import org.metagene.genestrip.util.StreamProvider;

public class KMerTrieFileGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;

	@SafeVarargs
	public KMerTrieFileGoal(GSProject project, String name, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			Goal<GSProject>... deps) {
		super(project, name, project.getTrieFile(), deps);
		this.taxNodesGoal = taxNodesGoal;
	}
	
	@Override
	protected void makeFile(File trieFile) {
		try {
			KMerTrie<String> trie = new KMerTrie<String>(1, getProject().getkMserSize(), false);

			Set<String> taxIds = new HashSet<String>();
			for (TaxIdNode node : taxNodesGoal.get()) {
				taxIds.add(node.getTaxId());
			}

			KrakenResultFastqMergeListener filter = KrakenResultFastqMergeListener.createFilterByTaxId(taxIds,
					KrakenResultFastqMergeListener.fillKMerTrie(trie, null));
			KrakenResultFastqMerger krakenKMerFastqMerger = new KrakenResultFastqMerger(
					getProject().getConfig().getMaxReadSizeBytes());

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Reading file " + getProject().getKrakenOutFile());
				getLogger().info("Reading file " + getProject().getKmerFastqFile());
			}
			
			InputStream stream1 = StreamProvider.getInputStreamForFile(getProject().getKrakenOutFile());
			InputStream stream2 = StreamProvider.getInputStreamForFile(getProject().getKmerFastqFile());			
			krakenKMerFastqMerger.process(stream1, stream1, filter);
			stream1.close();
			stream2.close();

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Trie entries: " + trie.getEntries());
				getLogger().info("File save " + trieFile);
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
}
