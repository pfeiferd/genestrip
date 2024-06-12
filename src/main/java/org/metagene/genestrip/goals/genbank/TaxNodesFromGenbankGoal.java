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
package org.metagene.genestrip.goals.genbank;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSConfigKey.SeqType;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.goals.refseq.RefSeqFnaFilesDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AbstractRefSeqFastaReader;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.StringLongDigitTrie;
import org.metagene.genestrip.util.StringLongDigitTrie.StringLong;

public class TaxNodesFromGenbankGoal extends ObjectGoal<Set<TaxIdNode>, GSProject> {
	private final ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;
	private final RefSeqFnaFilesDownloadGoal fnaFilesGoal;
	private final ObjectGoal<AccessionMap, GSProject> accessionMapGoal;

	@SafeVarargs
	public TaxNodesFromGenbankGoal(GSProject project, 
			ObjectGoal<Set<RefSeqCategory>[], GSProject> categoriesGoal,
			ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal, RefSeqFnaFilesDownloadGoal fnaFilesGoal,
			ObjectGoal<AccessionMap, GSProject> accessionMapGoal, Goal<GSProject>... deps) {
		super(project, GSGoalKey.TAXFROMGENBANK, Goal.append(deps, categoriesGoal, taxNodesGoal, fnaFilesGoal, accessionMapGoal));
		this.categoriesGoal = categoriesGoal;
		this.taxNodesGoal = taxNodesGoal;
		this.fnaFilesGoal = fnaFilesGoal;
		this.accessionMapGoal = accessionMapGoal;
	}

	@Override
	public void makeThis() {
		try {
			Set<TaxIdNode> missingTaxIds = new HashSet<TaxIdNode>();
			// We only get Genomic data from genbank (so far) - so if just RNA is wanted, there is no need to access it.
			if (!SeqType.RNA.equals(configValue(GSConfigKey.SEQ_TYPE))) {
				AbstractRefSeqFastaReader fastaReader = new AbstractRefSeqFastaReader(
						intConfigValue(GSConfigKey.FASTA_LINE_SIZE_BYTES), taxNodesGoal.get(), accessionMapGoal.get(),
						intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID)) {
					@Override
					protected void dataLine() throws IOException {
					}
				};

				for (File fnaFile : fnaFilesGoal.getFiles()) {
					RefSeqCategory cat = fnaFilesGoal.getCategoryForFile(fnaFile);
					if (categoriesGoal.get()[0].contains(cat)) {
						fastaReader.readFasta(fnaFile);
					}
				}
				StringLongDigitTrie trie = fastaReader.getRegionsPerTaxid();
				int limit = intConfigValue(GSConfigKey.REQ_SEQ_LIMIT_FOR_GENBANK);
				for (TaxIdNode node : taxNodesGoal.get()) {
					StringLong value = trie.get(node.getTaxId());
					long regions = value == null ? 0 : value.getLongValue();
					if (regions < limit) {
						missingTaxIds.add(node);
					}
				}
			}
			set(missingTaxIds);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}