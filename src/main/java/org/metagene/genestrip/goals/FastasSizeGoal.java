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
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.zip.GZIPInputStream;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.ObjectGoal;

public class FastasSizeGoal extends ObjectGoal<Long, GSProject> {
	private final FastaFileDownloadGoal fastaFileDownloadGoal;
	private final int maxCheckForGuess;

	public FastasSizeGoal(GSProject project, FastaFileDownloadGoal fastaDownloadGoal) {
		super(project, "fastassize");
		this.fastaFileDownloadGoal = fastaDownloadGoal;
		maxCheckForGuess = 10;
	}

	@Override
	public void makeThis() {
		int i = 0;
		long cSizeSum = 0;
		long uSizeSum = 0;

		for (File file : fastaFileDownloadGoal.getAvailableFiles()) {
			if (i == maxCheckForGuess) {
				break;
			}
			cSizeSum += file.length();
			uSizeSum += getUncompressedSize(file);
			i++;
		}
		if (i < maxCheckForGuess) {
			set(uSizeSum);
		} else {
			double cRate = ((double) uSizeSum) / cSizeSum;

			long sizeSum = 0;
			for (File file : fastaFileDownloadGoal.getAvailableFiles()) {
				sizeSum += file.length();
			}
			set((long) (sizeSum * cRate));
		}
	}

	protected long getUncompressedSize(File file) {
		try {
			ReadableByteChannel fc;
			fc = Channels.newChannel(new FileInputStream(file));
			GZIPInputStream gis = new GZIPInputStream(Channels.newInputStream(fc));
			long uSize;
			for (uSize = 0; gis.read() != -1; uSize++) {
			}
			return uSize;
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
