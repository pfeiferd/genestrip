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
import java.io.InputStream;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;

public class FastasSizeGoal extends ObjectGoal<Long, GSProject> {
	private final FastaFileDownloadGoal fastaFileDownloadGoal;
	private final int maxCheckForGuess;

	@SafeVarargs
	public FastasSizeGoal(GSProject project, String name, int maxCheckForGuess, FastaFileDownloadGoal fastaDownloadGoal,
			Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, fastaDownloadGoal));
		this.fastaFileDownloadGoal = fastaDownloadGoal;
		this.maxCheckForGuess = maxCheckForGuess;
	}

	@Override
	public void makeThis() {
		int i = 0;
		long cSizeSum = 0;
		long uSizeSum = 0;

		if (getLogger().isInfoEnabled()) {
			getLogger().info("Guessing total size via " + maxCheckForGuess + " fast.gz files.");
		}
		for (File file : fastaFileDownloadGoal.getAvailableFiles()) {
			if (i == maxCheckForGuess) {
				break;
			}
			cSizeSum += file.length();
			uSizeSum += getUncompressedSize(file);
			i++;
		}
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Total compressed size " + cSizeSum + " bytes, total uncompressed size " + uSizeSum + " bytes.");
		}
		if (i < maxCheckForGuess) {
			set(uSizeSum);
		} else {
			double cRate = ((double) uSizeSum) / cSizeSum;
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Estimated compression ratio is: " + cRate);
			}
			

			long sizeSum = 0;
			for (File file : fastaFileDownloadGoal.getAvailableFiles()) {
				sizeSum += file.length();
			}
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Total compressed size: " + (sizeSum / 1024 / 1024) + " MBytes");
			}			
			set((long) (sizeSum * cRate));
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Estimated total uncompressed size: " + (get() / 1024 / 1024) + " MBytes");
			}			
		}
	}

	protected long getUncompressedSize(File file) {
		try {
			InputStream stream = StreamProvider.getInputStreamForFile(file);
			long uSize;
			for (uSize = 0; stream.read() != -1; uSize++)
				;
			return uSize;
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
