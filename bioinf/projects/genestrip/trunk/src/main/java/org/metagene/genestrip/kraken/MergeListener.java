package org.metagene.genestrip.kraken;

import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;

import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.trie.KMerTrie;

public interface MergeListener {
	MergeListener STD_ERR_DEFAULT_UPDATE = new MergeListener() {
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

	public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read, byte[] readProbs);

	public static MergeListener fillKMerTrie(KMerTrie<String> trie, MergeListener delegate) {
		return new MergeListener() {
			@Override
			public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read,
					byte[] readProbs) {
				trie.put(read, 0, taxid);
				if (readCount % 10000 == 0) {
					if (KrakenKMerFastqMerger.logger.isInfoEnabled()) {
						KrakenKMerFastqMerger.logger.info("Trie entries:" + trie.getEntries());
						KrakenKMerFastqMerger.logger.info("Trie put ratio:" + ((double) trie.getEntries() / readCount));
					}
				}
				if (delegate != null) {
					delegate.newTaxidForRead(readCount, taxid, readDescriptor, read, readProbs);
				}
			}
		};
	}

	public static MergeListener createFilterByTaxIdNodes(Set<TaxIdNode> taxIdNodes, MergeListener delegate) {
		Set<String> taxIds = new HashSet<String>();
		for (TaxIdNode node : taxIdNodes) {
			taxIds.add(node.getTaxId());
		}
		return createExcludeTaxIds(taxIds, delegate);
	}

	public static MergeListener createFastQOutputFilterByTaxId(PrintStream printStream, MergeListener delegate) {
		return new MergeListener() {
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

	public static MergeListener createExcludeTaxIds(final Set<String> taxIds, final MergeListener delegate) {
		return new MergeListener() {
			@Override
			public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read,
					byte[] readProbs) {
				if (!taxIds.contains(taxid)) {
					delegate.newTaxidForRead(readCount, taxid, readDescriptor, read, readProbs);
				}
			}
		};
	}

	public static MergeListener createFilterByTaxId(final Set<String> taxIds, final MergeListener delegate) {
		return new MergeListener() {
			@Override
			public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read,
					byte[] readProbs) {
				if (taxIds.contains(taxid)) {
					delegate.newTaxidForRead(readCount, taxid, readDescriptor, read, readProbs);
				}
			}
		};
	}
}