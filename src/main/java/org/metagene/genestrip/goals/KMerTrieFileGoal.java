package org.metagene.genestrip.goals;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.fastqgen.KrakenKMerFastqMerger;
import org.metagene.genestrip.fastqgen.KrakenKMerFastqMerger.FilterListener;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.KMerTrie;

public class KMerTrieFileGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;

	@SafeVarargs
	public KMerTrieFileGoal(GSProject project, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, Goal<GSProject>... deps) {
		super(project, "triegen", project.getTrieFile(), deps);
		this.taxNodesGoal = taxNodesGoal;
	}

	@Override
	protected void makeFile(File trieFile) {
		try {
			KMerTrie<String> trie = new KMerTrie<String>(getProject().getkMserSize());

			Set<String> taxIds = new HashSet<String>();
			for (TaxIdNode node : taxNodesGoal.get()) {
				taxIds.add(node.getTaxId());
			}

			FilterListener filter = KrakenKMerFastqMerger.createFilterByTaxId(taxIds,
					KrakenKMerFastqMerger.fillKMerTrie(trie, null));
			KrakenKMerFastqMerger krakenKMerFastqMerger = new KrakenKMerFastqMerger(
					getProject().getConfig().getMaxReadSizeBytes());

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Reading file " + getProject().getKrakenOutFile());
				getLogger().info("Reading file " + getProject().getKmerFastqFile());
			}
			FileInputStream fStream = new FileInputStream(getProject().getKmerFastqFile());
			GZIPInputStream gStream = new GZIPInputStream(fStream, 4096);
			krakenKMerFastqMerger.process(new BufferedInputStream(new FileInputStream(getProject().getKrakenOutFile())),
					new BufferedInputStream(gStream), filter);

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
