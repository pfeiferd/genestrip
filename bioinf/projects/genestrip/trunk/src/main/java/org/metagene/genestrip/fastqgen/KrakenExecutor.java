package org.metagene.genestrip.fastqgen;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.text.MessageFormat;

public class KrakenExecutor {
	private final String binFolder;
	private final String execCommand;
	
	public KrakenExecutor(String binFolder, String execCommand) {
		this.execCommand = execCommand;
		this.binFolder = binFolder;
	}
	
	public String genExecLine(String database, File fastq) {
		return MessageFormat.format(execCommand, binFolder, database, fastq);
	}
	
	public void execute(String database, File fastq, File outFile) throws InterruptedException, IOException {
		Process process = Runtime.getRuntime().exec(genExecLine(database, fastq.getAbsoluteFile()));
		// TODO: Better write to gzipped stream?
		ReadableByteChannel readableByteChannel = Channels.newChannel(process.getInputStream());
		FileOutputStream out = new FileOutputStream(outFile);
		out.getChannel().transferFrom(readableByteChannel, 0, Long.MAX_VALUE);
		out.close();
		int res = process.waitFor();
		if (res != 0) {
			throw new IllegalStateException("Kraken terminated unsuccesfully");
		}
	}
}
