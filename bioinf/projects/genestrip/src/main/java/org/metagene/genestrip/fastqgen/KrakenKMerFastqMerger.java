package org.metagene.genestrip.fastqgen;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.KMerTrie;

public class KrakenKMerFastqMerger {
	public static FilterListener STD_ERR_DEFAULT_UPDATE = new FilterListener() {
		@Override
		public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read,
				byte[] readProbs) {
			System.err.print("Read: ");
			System.err.println(readCount);
			System.err.print("New taxid returned: ");
			System.err.println(taxid);
			System.err.println(readDescriptor);
			System.err.println(read);
			System.err.println(readProbs);
		}
	};

	public Map<String, Integer> process(InputStream bufferedInFromKraken, InputStream bufferedInFastQ,
			FilterListener update) throws IOException {
		byte[] krakenChars = new byte[2048];
		byte[] readDescriptor = new byte[2048];
		byte[] read = new byte[2 * 4096];
		byte[] readProbs = new byte[2 * 4096];

		DigitTrieNode root = new DigitTrieNode();

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

		Map<String, Integer> map = new HashMap<String, Integer>();
		root.collect(map);
		return map;
	}

	public static void print(Map<String, Integer> map, OutputStream out) {
		PrintStream pOut = new PrintStream(out);

		for (String taxid : map.keySet()) {
			pOut.print(taxid);
			pOut.print(';');
			pOut.println(map.get(taxid));
		}
	}

	// This Trie avoids generating useless key string objects when analyzing the
	// input stream.
	// It is very memory and time efficient...
	private class DigitTrieNode {
		private DigitTrieNode[] children;
		int value;
		String taxId;

		public DigitTrieNode() {
		}

		public String inc(byte[] seq, int start, int end) {
			int index;
			DigitTrieNode node = this, child;
			for (int i = start; i <= end; i++, node = child) {
				index = seq[i] - '0';
				if (node.children == null) {
					node.children = new DigitTrieNode[10];
				}
				child = node.children[index];
				if (child == null) {
					child = new DigitTrieNode();
					node.children[index] = child;
				}
			}
			node.value++;
			if (node.value == 1) {
				node.taxId = new String(seq, start, end - start + 1);
			}
			return node.taxId;
		}

		public void collect(Map<String, Integer> map) {
			if (value > 0) {
				map.put(taxId, value);
			}
			if (children != null) {
				for (int k = 0; k < children.length; k++) {
					if (children[k] != null) {
						children[k].collect(map);
					}
				}
			}
		}
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