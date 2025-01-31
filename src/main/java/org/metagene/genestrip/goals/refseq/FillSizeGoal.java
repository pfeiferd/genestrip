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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.StringLongDigitTrie.StringLong;

public class FillSizeGoal extends ObjectGoal<Long, GSProject> {
	private final ObjectGoal<Set<RefSeqCategory>, GSProject> categoriesGoal;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionMapGoal;

	@SafeVarargs
	public FillSizeGoal(GSProject project, ObjectGoal<Set<RefSeqCategory>, GSProject> categoriesGoal,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
			ObjectGoal<Map<File, TaxIdNode>, GSProject> additionalGoal,
			ObjectGoal<AccessionMap, GSProject> accessionMapGoal, Goal<GSProject>... deps) {
		super(project, GSGoalKey.FILLSIZE,
				Goal.append(deps, categoriesGoal, taxNodesGoal, fnaFilesGoal, accessionMapGoal, additionalGoal));
		this.categoriesGoal = categoriesGoal;
		this.taxNodesGoal = taxNodesGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.additionalGoal = additionalGoal;
		this.accessionMapGoal = accessionMapGoal;
	}

	@Override
	protected void doMakeThis() {
		try {
			MyFastaReader fastaReader = new MyFastaReader(intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES),
					taxNodesGoal.get(), accessionMapGoal.get(), intConfigValue(GSConfigKey.KMER_SIZE),
					intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID),
					(Rank) configValue(GSConfigKey.MAX_GENOMES_PER_TAXID_RANK),
					longConfigValue(GSConfigKey.MAX_KMERS_PER_TAXID));

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
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Store size determined in kmers: " + fastaReader.getCounter());
			}

			set(fastaReader.getCounter());
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected static class MyFastaReader extends AbstractRefSeqFastaReader {
		private long counter;

		public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k,
				int maxGenomesPerTaxId, Rank maxGenomesPerTaxIdRank, long maxKmersPerTaxId) {
			super(bufferSize, taxNodes, accessionMap, k, maxGenomesPerTaxId, maxGenomesPerTaxIdRank, maxKmersPerTaxId);
			counter = 0;
		}

		@Override
		protected void endRegion() {
			super.endRegion();
			if (includeRegion) {
				counter -= k - 1;
			}
		}

		@Override
		protected void dataLine() {
			super.dataLine();
			if (includeRegion) {
				counter += size - 1;
			}
		}

		public long getCounter() {
			return counter;
		}

		@Override
		protected void done() {
			if (getLogger().isTraceEnabled()) {
				List<StringLong> values = new ArrayList<StringLong>();
				regionsPerTaxid.collect(values);
				getLogger().trace("Included regions per taxid: " + values);
			}
			super.done();
		}
	}
}