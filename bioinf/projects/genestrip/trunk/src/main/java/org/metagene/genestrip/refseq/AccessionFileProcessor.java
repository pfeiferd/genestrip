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
package org.metagene.genestrip.refseq;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.util.BufferedLineReader;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.StreamProvider;
import org.metagene.genestrip.util.StreamProvider.ByteCountingInputStreamAccess;

public abstract class AccessionFileProcessor {
	private static int MAX_LINE_SIZE = 2048;

	protected static final String[] ALL_GENOMIC_ACCESSION_PREFIXES = { "AC_", "NC_", "NG_", "NT_", "NW_", "NZ_" };

	protected static final String[] COMPLETE_GENOMIC_ACCESSION_PREFIXES = { "AC_", "NC_", "NZ_" };

	protected final Log logger = LogFactory.getLog("accreader");

	private final long recordLogCycle = 1000 * 1000 * 10;
	private final boolean completeOnly;
	private final Collection<RefSeqCategory> categories;

	public AccessionFileProcessor(Collection<RefSeqCategory> categories, boolean completeOnly) {
		this.completeOnly = completeOnly;
		this.categories = categories;
	}

	public void processCatalog(File catalogFile) {
		try {
			ByteCountingInputStreamAccess byteCountAccess = StreamProvider
					.getByteCountingInputStreamForFile(catalogFile, false);

			long totalCatSize = Files.size(catalogFile.toPath());
			byte[] target = new byte[MAX_LINE_SIZE];

			// The file is huge, apache csv reader would be too slow and burn too many
			// string.
			// Therefore manual coding for parsing and
			int size;
			BufferedLineReader reader = new BufferedLineReader(byteCountAccess.getInputStream());

			long recordCounter = 0;
			long startTime = System.currentTimeMillis();
			while ((size = reader.nextLine(target)) > 0) {
				int pos1 = ByteArrayUtil.indexOf(target, 0, size, '\t');
				int pos2 = ByteArrayUtil.indexOf(target, pos1 + 1, size, '\t');
				int pos3 = ByteArrayUtil.indexOf(target, pos2 + 1, size, '\t');
				int pos4 = ByteArrayUtil.indexOf(target, pos3 + 1, size, '\t');
				if (isGenomicAccession(target, pos2 + 1)) {
					if (containsCategory(target, pos3 + 1, pos4, categories)) {
						handleEntry(target, pos1, pos2 + 1, pos3);
					}
				}

				if (recordCounter % recordLogCycle == 0) {
					if (logger.isInfoEnabled()) {
						double ratio = byteCountAccess.getBytesRead() / (double) totalCatSize;
						long stopTime = System.currentTimeMillis();

						double diff = (stopTime - startTime);
						double totalTime = diff / ratio;
						double totalHours = totalTime / 1000 / 60 / 60;

						logger.info("Elapsed hours:" + diff / 1000 / 60 / 60);
						logger.info("Estimated total hours:" + totalHours);
						logger.info("Records processed: " + recordCounter);
					}
				}
				recordCounter++;
			}
			reader.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected abstract void handleEntry(byte[] target, int taxIdEnd, int accessionStart, int accessionEnd);

	protected boolean containsCategory(byte[] outerArray, int start, int end, Collection<RefSeqCategory> categories) {
		for (RefSeqCategory cat : categories) {
			if (ByteArrayUtil.indexOf(outerArray, start, end, cat.getDirectory()) != -1) {
				return true;
			}
		}
		return false;
	}

	protected boolean isGenomicAccession(byte[] outerArray, int start) {
		for (String prefix : completeOnly ? COMPLETE_GENOMIC_ACCESSION_PREFIXES : ALL_GENOMIC_ACCESSION_PREFIXES) {
			if (ByteArrayUtil.startsWith(outerArray, start, prefix)) {
				return true;
			}
		}
		return false;
	}
}
