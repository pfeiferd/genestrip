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
	private final boolean writedFiltered;

	@SafeVarargs
	public ClassifyGoal(GSProject project, String name, File fastq, KMerTrieFileGoal trieGoal, boolean writeFiltered,
			Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, fastq, FileType.CSV, false), ArraysUtil.append(deps, trieGoal));
		this.fastq = fastq;
		this.trieGoal = trieGoal;
		this.writedFiltered = writeFiltered;
	}
	
	@Override
	protected void makeFile(File file) {
		try {
			File filteredFile = null;
			if (writedFiltered) {
				filteredFile = getProject().getOutputFile(getName(), fastq, FileType.FASTQ,false);
			}
			
			@SuppressWarnings("unchecked")
			KMerTrie<String> trie = (KMerTrie<String>) KMerTrie.load(trieGoal.getFile());
			Map<String, Long> res = new FastqTrieClassifier(trie, getProject().getConfig().getMaxReadSizeBytes())
					.runClassifier(fastq, filteredFile);
			PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
			CountingDigitTrie.print(res, out);
			out.close();
		} catch (IOException | ClassNotFoundException e) {
			throw new RuntimeException(e);
		}
	}

}
