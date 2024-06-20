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
package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractStoreFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class FillBloomFilterGoal extends ObjectGoal<MurmurCGATBloomFilter, GSProject> {
	private final ObjectGoal<Set<RefSeqCategory>, GSProject> categoriesGoal;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionMapGoal;
	private final ObjectGoal<Long, GSProject> sizeGoal;

	@SafeVarargs
	public FillBloomFilterGoal(GSProject project, ObjectGoal<Set<RefSeqCategory>, GSProject> categoriesGoal,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
			ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalGoal,
			ObjectGoal<AccessionMap, GSProject> accessionMapGoal, FillSizeGoal sizeGoal, Goal<GSProject>... deps) {
		super(project, GSGoalKey.TEMPINDEX, Goal.append(deps, taxNodesGoal, fnaFilesGoal, accessionMapGoal, sizeGoal, additionalGoal));
		this.categoriesGoal = categoriesGoal;
		this.taxNodesGoal = taxNodesGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.additionalGoal = additionalGoal;
		this.accessionMapGoal = accessionMapGoal;
		this.sizeGoal = sizeGoal;
	}
	
	@Override
	protected void allDependentsMade() {
		// To save memory...
		doCleanThis();
	}

	@Override
	protected void doMakeThis() {
		try {
			MurmurCGATBloomFilter filter = new MurmurCGATBloomFilter(intConfigValue(GSConfigKey.KMER_SIZE),
					doubleConfigValue(GSConfigKey.TEMP_BLOOM_FILTER_FPP));
			filter.ensureExpectedSize(sizeGoal.get(), false);

			MyFastaReader fastaReader = new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
					taxNodesGoal.get(), accessionMapGoal.get(), filter,
					intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID), intConfigValue(GSConfigKey.MAX_DUST));

			for (File fnaFile : fnaFilesGoal.getFiles()) {
				RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
				if (categoriesGoal.get().contains(cat)) {
					fastaReader.readFasta(fnaFile);
				}
			}
			Map<File, TaxIdNode> additionalMap = additionalGoal.get();
			for (File additionalFasta : additionalMap.keySet()) {
				TaxIdNode node = additionalMap.get(additionalFasta);
				if (taxNodesGoal.get().contains(node)) {
					fastaReader.ignoreAccessionMap(node);
					fastaReader.readFasta(additionalFasta);
				}
			}

			set(filter);
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Bloom filter entries: " + filter.getEntries());
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected static class MyFastaReader extends AbstractStoreFastaReader {
		private final MurmurCGATBloomFilter filter;

		public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap,
				MurmurCGATBloomFilter filter, int maxGenomesPerTaxId, int maxDust) {
			super(bufferSize, taxNodes, accessionMap, filter.getK(), maxGenomesPerTaxId, maxDust);
			this.filter = filter;
		}

		@Override
		protected void handleStore() {
			if (!filter.containsLong(byteRingBuffer.getKMer())) {
				filter.putLong(byteRingBuffer.getKMer());
			}
		}
	}
}