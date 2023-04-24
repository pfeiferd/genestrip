package org.metagene.genestrip.goals;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Map;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.trie.FastqTrieClassifier;
import org.metagene.genestrip.trie.KMerTrie;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.CountingDigitTrie;
import org.metagene.genestrip.util.StreamProvider;

public class ClassifyGoal extends FileListGoal<GSProject> {
	private final File fastq;
	private final KMerTrieFileGoal trieGoal;

	@SafeVarargs
	public ClassifyGoal(GSProject project, String name, File fastq, KMerTrieFileGoal trieGoal,
			Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, fastq, FileType.CSV, false), ArraysUtil.append(deps, trieGoal));
		this.fastq = fastq;
		this.trieGoal = trieGoal;
	}

	@Override
	protected void makeFile(File file) {
		try {
			@SuppressWarnings("unchecked")
			KMerTrie<String> trie = (KMerTrie<String>) KMerTrie.load(trieGoal.getOutputFile());
			Map<String, Long> res = new FastqTrieClassifier(trie, getProject().getConfig().getMaxReadSizeBytes())
					.runClassifier(fastq);
			PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
			CountingDigitTrie.print(res, out);
			out.close();
		} catch (IOException | ClassNotFoundException e) {
			throw new RuntimeException(e);
		}
	}

}
