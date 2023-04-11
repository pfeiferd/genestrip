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
import java.util.HashSet;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.KMerTrie;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.CountingDigitTrie;

public class KMerFastqTrieFileGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;

	@SafeVarargs
	public KMerFastqTrieFileGoal(GSProject project, String name, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			Goal<GSProject>... deps) {
		super(project, name, project.getTrieFile(), ArraysUtil.append(deps, taxNodesGoal));
		this.taxNodesGoal = taxNodesGoal;
	}

	@Override
	protected void makeFile(File trieFile) {
		try {
			CountingDigitTrie countingDigitTrie = new CountingDigitTrie();
			KMerTrie<String> trie = new KMerTrie<String>(1, getProject().getkMserSize(), false);

			Set<String> taxIds = new HashSet<String>();
			for (TaxIdNode node : taxNodesGoal.get()) {
				taxIds.add(node.getTaxId());
			}

			MyAbstractFastqReader fastqReader = new MyAbstractFastqReader(getProject().getConfig().getMaxReadSizeBytes()) {
				@Override
				protected void nextEntry() throws IOException {
					int pos;
					for (pos = readDescriptorSize - 1; pos >= 0; pos--) {
						if (readDescriptor[pos] == ':') {
							pos++;
							break;
						}
					}
					String taxid = countingDigitTrie.inc(readDescriptor, pos, readDescriptorSize - 1);
					if (taxIds.contains(taxid)) {
						trie.put(read, 0, taxid, false);
					}
				}
			};
			fastqReader.readFastq(getProject().getFilteredKmerFastqFile());

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Trie entries: " + trie.getEntries());
				getLogger().info("Saving file " + trieFile);
			}
			trie.compress();
			trie.save(trieFile);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Saved file " + trieFile);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	
	protected static abstract class MyAbstractFastqReader extends AbstractFastqReader {
		public MyAbstractFastqReader(int maxReadSizeBytes) {
			super(maxReadSizeBytes);
		}

		@Override
		public void readFastq(File file) throws IOException {
			super.readFastq(file);
		}
	}
}
