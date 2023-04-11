package org.metagene.genestrip.fastq;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.util.BufferedLineReader;
import org.metagene.genestrip.util.StreamProvider;

public abstract class AbstractFastqReader {
	protected final Log logger = LogFactory.getLog(getClass());

	protected final byte[] readDescriptor;
	protected int readDescriptorSize;

	protected final byte[] read;
	protected int readSize;

	protected final byte[] readProbs;
	protected int readProbsSize;

	protected long reads;

	private final BufferedLineReader bufferedLineReaderFastQ;

	public AbstractFastqReader(int maxReadSizeBytes) {
		readDescriptor = new byte[maxReadSizeBytes];
		read = new byte[maxReadSizeBytes];
		readProbs = new byte[maxReadSizeBytes];

		bufferedLineReaderFastQ = new BufferedLineReader();
	}

	protected Log getLogger() {
		return logger;
	}

	protected void readFastq(File file) throws IOException {
		if (logger.isInfoEnabled()) {
			logger.info("Reading fastq file " + file);
		}
		InputStream inputStream = StreamProvider.getInputStreamForFile(file);
		readFastq(inputStream);
		inputStream.close();
		if (logger.isInfoEnabled()) {
			logger.info("Closed fastq file " + file);
		}
	}

	protected void readFastq(InputStream inputStream) throws IOException {
		reads = 0;
		bufferedLineReaderFastQ.setInputStream(inputStream);

		start();
		for (readDescriptorSize = bufferedLineReaderFastQ
				.nextLine(readDescriptor); readDescriptorSize > 0; readDescriptorSize = bufferedLineReaderFastQ
						.nextLine(readDescriptor)) {
			readDescriptor[readDescriptorSize - 1] = 0;
			readSize = bufferedLineReaderFastQ.nextLine(read);
			read[readSize - 1] = 0;
			bufferedLineReaderFastQ.skipLine(); // Ignoring line 3.
			readProbsSize = bufferedLineReaderFastQ.nextLine(readProbs);
			readProbs[readProbsSize - 1] = 0;

			reads++;
			nextEntry();
		}
		done();
		if (logger.isInfoEnabled()) {
			logger.info("Total number of reads: " + reads);
		}
	}

	protected abstract void nextEntry() throws IOException;
	
	protected void done() throws IOException {};
	
	protected void start() throws IOException {};
}
