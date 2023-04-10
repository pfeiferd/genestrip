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

import java.io.File;
import java.io.IOException;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.fastqgen.KMerFastqGenerator;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.util.ArraysUtil;

public class KMerFastqGoal extends FileListGoal<GSProject> {
	private long addedKmers;
	private final FastasSizeGoal fastasSizeGoal;
	private final FastaFileDownloadGoal fastaDownloadGoal;

	@SafeVarargs
	public KMerFastqGoal(GSProject project, String name, FastasSizeGoal fastasSizeGoal,
			FastaFileDownloadGoal fastaDownloadGoal, Goal<GSProject>... deps) {
		super(project, name, project.getKmerFastqFile(), ArraysUtil.append(deps, fastasSizeGoal, fastaDownloadGoal));
		this.fastasSizeGoal = fastasSizeGoal;
		this.fastaDownloadGoal = fastaDownloadGoal;
	}

	@Override
	protected void makeFile(File fastq) {
		try {
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Generating fastq file " + fastq);
			}
			KMerFastqGenerator generator = new KMerFastqGenerator(getProject().getkMserSize(),
					getProject().getConfig().getMaxReadSizeBytes());
			addedKmers = generator.run(fastaDownloadGoal.getAvailableFiles(), fastq, getProject().getName(),
					fastasSizeGoal.get());
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Entered K-mers: " + addedKmers);
				getLogger().info("Generated fastq file " + fastq);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public long getAddedKmers() {
		return addedKmers;
	}
}
