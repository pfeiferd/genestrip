package org.metagene.genestrip.gen.goals;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.metagene.genestrip.fastqgen.KrakenKMerFastqMerger;
import org.metagene.genestrip.fastqgen.KrakenKMerFastqMerger.FilterListener;
import org.metagene.genestrip.gen.Project;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.KMerTrie;

public class KMerTrieFileGoal extends FileListGoal {
	private final Project project;
	private final ObjectGoal<Set<TaxIdNode>> taxNodesGoal;

	public KMerTrieFileGoal(Project project, ObjectGoal<Set<TaxIdNode>> taxNodesGoal, Goal... deps) {
		super("triegen", project.getTrieFile(), deps);
		this.project = project;
		this.taxNodesGoal = taxNodesGoal;
	}

	@Override
	protected void makeFile(File trieFile) {
		try {
			KMerTrie<String> trie = new KMerTrie<String>(project.getkMserSize());

			Set<String> taxIds = new HashSet<String>();
			for (TaxIdNode node : taxNodesGoal.get()) {
				taxIds.add(node.getTaxId());
			}

			FilterListener filter = KrakenKMerFastqMerger.createFilterByTaxId(taxIds, KrakenKMerFastqMerger.fillKMerTrie(trie, null));
			KrakenKMerFastqMerger krakenKMerFastqMerger = new KrakenKMerFastqMerger();

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Reading file " + project.getKrakenOutFile());
				getLogger().info("Reading file " + project.getKmerFastqFile());
			}
			FileInputStream fStream = new FileInputStream(project.getKmerFastqFile());
			GZIPInputStream gStream = new GZIPInputStream(fStream, 4096);
			krakenKMerFastqMerger.process(new BufferedInputStream(new FileInputStream(project.getKrakenOutFile())),
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
