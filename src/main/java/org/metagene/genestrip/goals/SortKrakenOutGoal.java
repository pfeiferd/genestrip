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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.kraken.SortExecutor;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;

public class SortKrakenOutGoal extends FileListGoal<GSProject> {
	private final Map<File, File> goalFileToKrakeFile;

	@SafeVarargs
	public SortKrakenOutGoal(GSProject project, String name, FileGoal<GSProject> krakenOutGoal,
			Goal<GSProject>... deps) {
		super(project, name, (List<File>)null, ArraysUtil.append(deps, krakenOutGoal));
		goalFileToKrakeFile = new HashMap<File, File>();
		for (File file : krakenOutGoal.getFiles()) {
			File goalFile = project.getOutputFile(name, file, FileType.KRAKEN_OUT);
			addFile(goalFile);
			goalFileToKrakeFile.put(goalFile, file);
		}
	}

	@Override
	protected void makeFile(File outFile) {
		SortExecutor sortExecutor = new SortExecutor(getProject().getConfig().getSortBin(),
				getProject().getConfig().getSortExecExpr());
		try {
			if (getLogger().isInfoEnabled()) {
				String execLine = sortExecutor.genExecLine(goalFileToKrakeFile.get(outFile));
				getLogger().info("Running sort with " + execLine);
			}
			OutputStream out = StreamProvider.getOutputStreamForFile(outFile);
			sortExecutor.execute(goalFileToKrakeFile.get(outFile), out);
			out.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		}
	}
}