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
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.goals.TrieFromKrakenResGoal.TaxidWithCount;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.KMerTrie;
import org.metagene.genestrip.store.KMerTrie.KMerTrieVisitor;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.StreamProvider;

public class KrakenResErrorGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<KMerTrie<TaxidWithCount>, GSProject> trieFromKrakenResGoal;

	@SafeVarargs
	public KrakenResErrorGoal(GSProject project, String name, File fastq,
			ObjectGoal<KMerTrie<TaxidWithCount>, GSProject> trieFromKrakenResGoal, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile("krakenerr", fastq, FileType.CSV, false),
				ArraysUtil.append(deps, trieFromKrakenResGoal));
		this.trieFromKrakenResGoal = trieFromKrakenResGoal;
	}

	@Override
	protected void makeFile(File file) throws IOException {
		final PrintStream ps = new PrintStream(StreamProvider.getOutputStreamForFile(file));

		ps.println("taxid;count;kmer");
		KMerTrie<TaxidWithCount> trie = trieFromKrakenResGoal.get();
		Map<String, Integer> errPerTaxid = new HashMap<String, Integer>();

		trie.visit(new KMerTrieVisitor<TrieFromKrakenResGoal.TaxidWithCount>() {
			@Override
			public void nextValue(KMerTrie<TaxidWithCount> trie, byte[] kmer, TaxidWithCount value) {
				if (value != null) {
					ps.print(value.getTaxid());
					ps.print(';');
					ps.print(value.getCount());
					ps.print(';');
					ByteArrayUtil.print(kmer, ps);
					ps.println(';');
					Integer res = errPerTaxid.get(value.getTaxid());
					res = (res == null) ? value.getCount() : res + value.getCount();
					errPerTaxid.put(value.getTaxid(), res);
				}
			}
		}, false);
		ps.println("taxid;aggregated count");
		for (Entry<String, Integer> entry : errPerTaxid.entrySet()) {
			ps.print(entry.getKey());
			ps.print(';');
			ps.print(entry.getValue());
			ps.println(';');
		}

		ps.close();
	}
}
