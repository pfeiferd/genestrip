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
package org.metagene.genestrip.goals;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.StreamProvider;

public class KrakenFastqFileGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;

	@SafeVarargs
	public KrakenFastqFileGoal(GSProject project, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			Goal<GSProject>... deps) {
		super(project, "krakenfastq", project.getFilteredKmerFastqFile(), deps);
		this.taxNodesGoal = taxNodesGoal;
	}

	@Override
	protected void makeFile(File fastq) {
		try {
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Writing file " + fastq);
			}
			PrintStream filteredOut = new PrintStream(StreamProvider.getOutputStreamForFile(fastq));

			KrakenResultFastqMergeListener filter = KrakenResultFastqMergeListener.createFilterByTaxIdNodes(
					taxNodesGoal.get(),
					KrakenResultFastqMergeListener.createFastQOutputFilterByTaxId(filteredOut, null));
			KrakenResultFastqMerger krakenKMerFastqMerger = new KrakenResultFastqMerger(
					getProject().getConfig().getMaxReadSizeBytes());

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Reading file " + getProject().getKrakenOutFile());
				getLogger().info("Reading file " + getProject().getKmerFastqFile());
			}
			krakenKMerFastqMerger.process(
					new BufferedInputStream(StreamProvider.getInputStreamForFile(getProject().getKrakenOutFile())),
					StreamProvider.getInputStreamForFile(getProject().getKmerFastqFile()), filter);

			filteredOut.close();
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Written file " + getProject().getFilteredKmerFastqFile());
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

}
