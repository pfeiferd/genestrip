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
import org.metagene.genestrip.trie.KMerTrie;
import org.metagene.genestrip.trie.KMerTrie.KMerTrieVisitor;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.StreamProvider;

public class KrakenResErrorGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<KMerTrie<TaxidWithCount>, GSProject> trieFromKrakenResGoal;

	@SafeVarargs
	public KrakenResErrorGoal(GSProject project, String name, File fastq,
			ObjectGoal<KMerTrie<TaxidWithCount>, GSProject> trieFromKrakenResGoal, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile("krakenerr", fastq, FileType.CSV),
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
