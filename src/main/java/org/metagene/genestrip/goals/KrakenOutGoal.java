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
import java.io.OutputStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.kraken.KrakenExecutor;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.util.StreamProvider;

public class KrakenOutGoal extends FileListGoal<GSProject> {
	private final Map<File, File> outFileTofastqFile;

	@SafeVarargs
	public KrakenOutGoal(GSProject project, String name, File fastqFile, Goal<GSProject>... deps) {
		this(project, name, Collections.singletonList(
				project.getOutputFile(name, fastqFile, FileType.KRAKEN_OUT, !project.getConfig().isUseKraken1())),
				deps);
	}

	@SafeVarargs
	public KrakenOutGoal(GSProject project, String name, List<File> fastqFiles, Goal<GSProject>... deps) {
		super(project, name, (List<File>) null, deps);

		outFileTofastqFile = new HashMap<File, File>();
		for (File fastqFile : fastqFiles) {
			File outFile = project.getOutputFile(name, fastqFile, FileType.KRAKEN_OUT,
					!project.getConfig().isUseKraken1());
			outFileTofastqFile.put(outFile, fastqFile);
			addFile(outFile);
		}
	}

	@Override
	protected void makeFile(File krakenOut) {
		KrakenExecutor krakenExecutor = new KrakenExecutor(getProject().getConfig().getKrakenBin(),
				getProject().getConfig().getKrakenExecExpr());
		try {
			File fastqFile = outFileTofastqFile.get(krakenOut);
			if (getLogger().isInfoEnabled()) {
				String execLine = krakenExecutor.genExecLine(getProject().getKrakenDB(), fastqFile, krakenOut);
				getLogger().info("Running kraken with " + execLine);
			}
			if (krakenExecutor.isWithFileForOutput()) {
				krakenExecutor.execute2(getProject().getKrakenDB(), fastqFile, krakenOut, System.out, System.err);
			}
			else {
				OutputStream out = StreamProvider.getOutputStreamForFile(krakenOut);
				krakenExecutor.execute2(getProject().getKrakenDB(), fastqFile, krakenOut, out, System.err);
				out.close();
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Finished kraken");
		}
	}
}
