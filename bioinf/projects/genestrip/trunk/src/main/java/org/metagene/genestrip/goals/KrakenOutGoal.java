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
import org.metagene.genestrip.kraken.KrakenExecutor;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;

public class KrakenOutGoal extends FileListGoal<GSProject> {
	private final File fastqFile;

	@SafeVarargs
	public KrakenOutGoal(GSProject project, String name, File fastqFile, File outfile,
			Goal<GSProject>... dependencies) {
		super(project, name, outfile, dependencies);
		this.fastqFile = fastqFile;
	}

	@Override
	protected void makeFile(File krakenOut) {

		KrakenExecutor krakenExecutor = new KrakenExecutor(getProject().getConfig().getKrakenBinFolder(),
				getProject().getConfig().getKrakenExecExpr());
		if (getLogger().isInfoEnabled()) {
			String execLine = krakenExecutor.genExecLine(getProject().getKrakenDB(), fastqFile);
			getLogger().info("Run kraken with " + execLine);
		}
		try {
			krakenExecutor.execute(getProject().getKrakenDB(), fastqFile, krakenOut);
		} catch (InterruptedException | IOException e) {
			throw new RuntimeException(e);
		}
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Finished kraken");
		}
	}
}
