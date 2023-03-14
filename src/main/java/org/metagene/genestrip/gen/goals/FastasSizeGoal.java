package org.metagene.genestrip.gen.goals;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.zip.GZIPInputStream;

import org.metagene.genestrip.make.ObjectGoal;

public class FastasSizeGoal extends ObjectGoal<Long> {
	private final FastaFileDownloadGoal fastaFileDownloadGoal;
	private final int maxCheckForGuess;

	public FastasSizeGoal(FastaFileDownloadGoal fastaDownloadGoal) {
		super("fastassize");
		this.fastaFileDownloadGoal = fastaDownloadGoal;
		maxCheckForGuess = 10;
	}

	@Override
	public void makeThis() {
		int i = 0;
		long cSizeSum = 0;
		long uSizeSum = 0;

		for (File file : fastaFileDownloadGoal.getFiles()) {
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
			for (File file : fastaFileDownloadGoal.getFiles()) {
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
