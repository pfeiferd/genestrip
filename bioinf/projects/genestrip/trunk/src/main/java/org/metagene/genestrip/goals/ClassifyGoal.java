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
import java.io.PrintStream;
import java.util.Map;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.trie.KMerTrie;
import org.metagene.genestrip.trie.FastqTrieClassifier;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.CountingDigitTrie;
import org.metagene.genestrip.util.StreamProvider;

public class ClassifyGoal extends FileListGoal<GSProject> {
	private final File fastq;
	private final KMerTrieFileGoal trieGoal;
	private final boolean writedFiltered;

	@SafeVarargs
	public ClassifyGoal(GSProject project, String name, File fastq, KMerTrieFileGoal trieGoal, boolean writeFiltered,
			Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, fastq, FileType.CSV, false),
				ArraysUtil.append(deps, trieGoal));
		this.fastq = fastq;
		this.trieGoal = trieGoal;
		this.writedFiltered = writeFiltered;
	}

	@Override
	protected void makeFile(File file) {
		FastqTrieClassifier c = null;
		try {
			File filteredFile = null;
			if (writedFiltered) {
				filteredFile = getProject().getOutputFile(getName(), fastq, FileType.FASTQ, false);
			}

			@SuppressWarnings("unchecked")
			KMerTrie<String> trie = (KMerTrie<String>) KMerTrie.load(trieGoal.getFile());

			c = new FastqTrieClassifier(trie, getProject().getConfig().getMaxReadSizeBytes(),
					getProject().getConfig().getThreadQueueSize(), getProject().getConfig().getThreads());
			Map<String, Long> res = c.runClassifier(fastq, filteredFile);
			c.dump();

			PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
			CountingDigitTrie.print(res, out);
			out.close();
		} catch (IOException | ClassNotFoundException e) {
			throw new RuntimeException(e);
		} finally {
			if (c != null) {
				c.dump();
			}
		}
	}
}
