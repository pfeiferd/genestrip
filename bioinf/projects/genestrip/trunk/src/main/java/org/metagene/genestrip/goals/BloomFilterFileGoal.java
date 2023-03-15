package org.metagene.genestrip.goals;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.bloom.KMerBloomIndex;
import org.metagene.genestrip.fastqgen.KrakenKMerFastqMerger;
import org.metagene.genestrip.fastqgen.KrakenKMerFastqMerger.FilterListener;
import org.metagene.genestrip.goals.KMerFastqGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class BloomFilterFileGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final KMerFastqGoal kMerFastqGoal;

	@SafeVarargs
	public BloomFilterFileGoal(GSProject project, KMerFastqGoal kMerFastqGoal,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, Goal<GSProject>... deps) {
		super(project, "bloomgen", project.getBloomFilterFile(), deps);
		this.taxNodesGoal = taxNodesGoal;
		this.kMerFastqGoal = kMerFastqGoal;
	}

	@Override
	protected void makeFile(File bloomFilterFile) {
		try {
			KMerBloomIndex bloomIndex = new KMerBloomIndex(bloomFilterFile.getName(), getProject().getkMserSize(),
					kMerFastqGoal.getAddedKmers(), 0.0001, null);

			Set<String> taxIds = new HashSet<String>();
			for (TaxIdNode node : taxNodesGoal.get()) {
				taxIds.add(node.getTaxId());
			}

			FilterListener filter = KrakenKMerFastqMerger.createFilterByTaxId(taxIds, new FilterListener() {
				@Override
				public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read,
						byte[] readProbs) {
					bloomIndex.putDirectKMer(read, 0);
				}
			});
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
				getLogger().info("File save " + bloomFilterFile);
			}
			bloomIndex.save(bloomFilterFile);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("File saved " + bloomFilterFile);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

}
