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
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.goals.AdditionalFastasGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class FillSizeGoal extends ObjectGoal<Long, GSProject> {
	private final ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final AdditionalFastasGoal additionalGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionMapGoal;

	@SafeVarargs
	public FillSizeGoal(GSProject project, String name, ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
			AdditionalFastasGoal additionalGoal, ObjectGoal<AccessionMap, GSProject> accessionMapGoal,
			Goal<GSProject>... deps) {
		super(project, name, Goal.append(deps, categoriesGoal, taxNodesGoal, fnaFilesGoal, accessionMapGoal));
		this.categoriesGoal = categoriesGoal;
		this.taxNodesGoal = taxNodesGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.additionalGoal = additionalGoal;
		this.accessionMapGoal = accessionMapGoal;
	}

	@Override
	public void makeThis() {
		try {
			Set<TaxIdNode> taxNodes = taxNodesGoal.get();
			MyFastaReader fastaReader = new MyFastaReader(getProject().getConfig().getMaxReadSizeBytes(),
					taxNodes, accessionMapGoal.get(), getProject().getConfig().getKMerSize());

			for (File fnaFile : fnaFilesGoal.getFiles()) {
				RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
				if (categoriesGoal.get()[0].contains(cat)) {
					fastaReader.readFasta(fnaFile);
				}
			}
			for (File additionalFasta : additionalGoal.getFiles()) {
				TaxIdNode additionalNode = additionalGoal.getTaxNodeForFile(additionalFasta);
				if (taxNodes.contains(additionalNode)) {
					fastaReader.ignoreAccessionMap(additionalNode);
					fastaReader.readFasta(additionalFasta);
				}
			}
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Store size determined: " + fastaReader.getCounter());
			}

			set(fastaReader.getCounter());
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected static class MyFastaReader extends AbstractRefSeqFastaReader {
		private long counter;
		private final int k;

		public MyFastaReader(int bufferSize, Set<TaxIdNode> taxNodes, AccessionMap accessionMap, int k) {
			super(bufferSize, taxNodes, accessionMap);
			this.k = k;
			counter = 0;
		}

		@Override
		protected void infoLine() throws IOException {
			if (includeRegion) {
				counter -= k - 1;
			}
			super.infoLine();
		}

		@Override
		protected void dataLine() throws IOException {
			if (includeRegion) {
				counter += size - 1;
			}
		}

		public long getCounter() {
			return counter;
		}

		@Override
		protected void done() throws IOException {
			if (includeRegion) {
				counter -= k - 1;
			}
		}
	}
}