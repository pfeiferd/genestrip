package org.metagene.genestrip.fastqgen;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.KMerTrie;
import org.metagene.genestrip.util.CountingDigitTrie;

public class KrakenKMerFastqMerger {
	public static FilterListener STD_ERR_DEFAULT_UPDATE=new FilterListener(){@Override public void newTaxidForRead(long readCount,String taxid,byte[]readDescriptor,byte[]read,byte[]readProbs){System.err.print("Read: ");System.err.println(readCount);System.err.print("New taxid returned: ");System.err.println(taxid);System.err.println(readDescriptor);System.err.println(read);System.err.println(readProbs);}};

	protected static final Log logger = LogFactory.getLog(KrakenKMerFastqMerger.class);

	private final byte[] krakenChars;
	private final byte[] readDescriptor;
	private final byte[] read;
	private final byte[] readProbs;

	public KrakenKMerFastqMerger(int maxReadSizeBytes) {
		krakenChars = new byte[2048];
		readDescriptor = new byte[2048];
		read = new byte[maxReadSizeBytes];
		readProbs = new byte[maxReadSizeBytes];
	}

	public Map<String, Long> process(InputStream bufferedInFromKraken, InputStream bufferedInFastQ,
			FilterListener update) throws IOException {
		CountingDigitTrie root = new CountingDigitTrie();

		int krakenPos;
		long readCount = 0;
		for (int c = bufferedInFromKraken.read(), d = bufferedInFastQ.read(); c != -1
				&& d != -1; c = bufferedInFromKraken.read(), d = bufferedInFastQ.read()) {
			for (krakenPos = 0; c != -1 && c != '\n'; krakenPos++) {
				krakenChars[krakenPos] = (byte) c;
				c = bufferedInFromKraken.read();
			}

			for (int pos = 0; d != -1 && d != '\n'; pos++) {
				readDescriptor[pos] = (byte) d;
				d = bufferedInFastQ.read();
			}

			d = bufferedInFastQ.read();
			for (int pos = 0; d != -1 && d != '\n'; pos++) {
				read[pos] = (byte) d;
				d = bufferedInFastQ.read();
			}
			// Ignore line with "+...":
			d = bufferedInFastQ.read();
			for (; d != -1 && d != '\n';) {
				d = bufferedInFastQ.read();
			}
			d = bufferedInFastQ.read();
			for (int pos = 0; d != -1 && d != '\n'; pos++) {
				readProbs[pos] = (byte) d;
				d = bufferedInFastQ.read();
			}
			readCount++;
			krakenPos--;
			int start;
			int end = krakenPos - 2;
			for (start = krakenPos; start > 0; start--) {
				if (krakenChars[start] == ':') {
					end = krakenPos;
				}
				if (krakenChars[start] == '\t') {
					break;
				}
			}
			String taxid = root.inc(krakenChars, start + 1, end - 2);

			int i;
			for (i = 1; i < krakenChars.length && krakenChars[i + 1] != '\t'; i++) {
				if (krakenChars[i + 1] != readDescriptor[i]) {
					throw new IllegalStateException("In consistent files for read " + readCount);
				}
			}
			if (i == krakenChars.length) {
				throw new IllegalStateException("In consistent kraken output...");
			}

			if (taxid != null && update != null) {
				update.newTaxidForRead(readCount, taxid, readDescriptor, read, readProbs);
			}
		}
		bufferedInFastQ.close();
		bufferedInFromKraken.close();

		Map<String, Long> map = new HashMap<String, Long>();
		root.collect(map);
		return map;
	}

	public static FilterListener createFilterByTaxIdNodes(Set<TaxIdNode> taxIdNodes, FilterListener delegate) {
		Set<String> taxIds = new HashSet<String>();
		for (TaxIdNode node : taxIdNodes) {
			taxIds.add(node.getTaxId());
		}
		return createExcludeTaxIds(taxIds, delegate);
	}

	public static FilterListener createFilterByTaxId(final Set<String> taxIds, final FilterListener delegate) {
		return new FilterListener() {
			@Override
			public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read,
					byte[] readProbs) {
				if (taxIds.contains(taxid)) {
					delegate.newTaxidForRead(readCount, taxid, readDescriptor, read, readProbs);
				}
			}
		};
	}

	public static FilterListener createExcludeTaxIds(final Set<String> taxIds, final FilterListener delegate) {
		return new FilterListener() {
			@Override
			public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read,
					byte[] readProbs) {
				if (!taxIds.contains(taxid)) {
					delegate.newTaxidForRead(readCount, taxid, readDescriptor, read, readProbs);
				}
			}
		};
	}

	public static FilterListener createFastQOutputFilterByTaxId(PrintStream printStream, FilterListener delegate) {
		return new FilterListener() {
			@Override
			public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read,
					byte[] readProbs) {
				printStream.print(readDescriptor);
				printStream.print(":taxid:");
				printStream.println(taxid);
				printStream.println(read);
				printStream.println("+");
				printStream.println(readProbs);

				if (delegate != null) {
					delegate.newTaxidForRead(readCount, taxid, readDescriptor, read, readProbs);
				}
			}
		};
	}

	public static FilterListener fillKMerTrie(KMerTrie<String> trie, FilterListener delegate) {
		return new FilterListener() {
			@Override
			public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read,
					byte[] readProbs) {
				trie.put(read, 0, taxid);
				if (readCount % 10000 == 0) {
					if (logger.isInfoEnabled()) {
						logger.info("Trie entries:" + trie.getEntries());
						logger.info("Trie put ratio:" + ((double) trie.getEntries() / readCount));
					}
				}
				if (delegate != null) {
					delegate.newTaxidForRead(readCount, taxid, readDescriptor, read, readProbs);
				}
			}
		};
	}

	public interface FilterListener {
		public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read, byte[] readProbs);
	}
}
