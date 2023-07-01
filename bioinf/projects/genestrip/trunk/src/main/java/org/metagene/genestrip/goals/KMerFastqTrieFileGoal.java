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
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.fastq.AbstractFastqReader;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.KMerStore;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.CountingDigitTrie;

public class KMerFastqTrieFileGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final KrakenFastqFileGoal krakenFastqGoal;

	@SafeVarargs
	public KMerFastqTrieFileGoal(GSProject project, String name, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			KrakenFastqFileGoal krakenFastqGoal, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.SER),
				ArraysUtil.append(deps, taxNodesGoal, krakenFastqGoal));
		this.taxNodesGoal = taxNodesGoal;
		this.krakenFastqGoal = krakenFastqGoal;
	}

	@Override
	protected void makeFile(File trieFile) {
		try {
			CountingDigitTrie countingDigitTrie = new CountingDigitTrie();
			KMerStore<String> trie = getProject().getKMerStoreFactory().createKMerStore(String.class);

			Set<String> taxIds = new HashSet<String>();
			for (TaxIdNode node : taxNodesGoal.get()) {
				taxIds.add(node.getTaxId());
			}

			MyAbstractFastqReader fastqReader = new MyAbstractFastqReader(getProject().getConfig().getkMerSize(),
					getProject().getConfig().getMaxReadSizeBytes()) {
				@Override
				protected void nextEntry(ReadEntry readEntry) throws IOException {
					int pos;
					for (pos = readEntry.readDescriptorSize - 1; pos >= 0; pos--) {
						if (readEntry.readDescriptor[pos] == ':') {
							pos++;
							break;
						}
					}
					String taxid = countingDigitTrie.inc(readEntry.readDescriptor, pos, readEntry.readDescriptorSize);
					if (taxIds.contains(taxid)) {
						if (!trie.put(readEntry.read, 0, taxid, false)) {
							if (getLogger().isWarnEnabled()) {
								getLogger().warn("Duplicate entry for read regarding taxid " + taxid);
							}							
						}
					}
				}
			};
			for (File file : krakenFastqGoal.getFiles()) {
				fastqReader.readFastq(file);
			}

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Trie entries: " + trie.getEntries());
				getLogger().info("Saving file " + trieFile);
			}
			trie.optimize();
			KMerStore.save(trie, trieFile);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Saved file " + trieFile);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected static abstract class MyAbstractFastqReader extends AbstractFastqReader {
		public MyAbstractFastqReader(int k, int maxReadSizeBytes) {
			super(k, maxReadSizeBytes, 0, 0);
		}

		@Override
		public void readFastq(File file) throws IOException {
			super.readFastq(file);
		}
	}
}
