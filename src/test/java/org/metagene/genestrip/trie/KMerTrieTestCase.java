package org.metagene.genestrip.trie;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Random;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.Main;
import org.metagene.genestrip.fastqgen.KrakenKMerFastqMerger;
import org.metagene.genestrip.fastqgen.KrakenKMerFastqMerger.FilterListener;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.CountingDigitTrie;

import com.google.common.hash.BloomFilter;
import com.google.common.hash.Funnel;
import com.google.common.hash.PrimitiveSink;

import junit.framework.TestCase;

public class KMerTrieTestCase extends TestCase {
	public void testPutGet() throws IOException {
		Main main = new Main();
		main.parseAndRun(new String[] { "bart_h" });

		GSProject project = main.getProject();

		File bartHReads = project.getKmerFastqFile();
		FileInputStream fStream = new FileInputStream(bartHReads);
		GZIPInputStream gStream = new GZIPInputStream(fStream, 4096);

		File fromKraken = project.getKrakenOutFile();

		KMerTrie<String> trie = new KMerTrie<String>(project.getkMserSize());

		@SuppressWarnings("serial")
		Funnel<byte[]> funnel = new Funnel<byte[]>() {
			@Override
			public void funnel(byte[] from, PrimitiveSink into) {
				for (int i = 0; i < trie.getLen(); i++) {
					into.putByte((byte) from[i]);
				}
			}
		};
		BloomFilter<byte[]> bloomFilter = BloomFilter.create(funnel, 5 * 1000 * 1000, 0.00001);

		@SuppressWarnings("unchecked")
		ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal = (ObjectGoal<Set<TaxIdNode>, GSProject>) main.getGenerator()
				.getGoal("taxids");
		Set<TaxIdNode> nodes = taxNodesGoal.get();

		FilterListener filter = KrakenKMerFastqMerger.createFilterByTaxIdNodes(nodes,
				KrakenKMerFastqMerger.fillKMerTrie(trie, new FilterListener() {
					@Override
					public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read,
							byte[] readProbs) {
						bloomFilter.put(read);
					}
				}));

		KrakenKMerFastqMerger krakenFilter = new KrakenKMerFastqMerger(project.getConfig().getMaxReadSizeBytes());

		CountingDigitTrie.print(krakenFilter.process(new BufferedInputStream(new FileInputStream(fromKraken)),
				new BufferedInputStream(gStream), filter), System.out);

		// Test uncompressed:
		checkTrie(fromKraken, bartHReads, krakenFilter, trie, bloomFilter, nodes);

		trie.compress();

		// Test compressed;
		checkTrie(fromKraken, bartHReads, krakenFilter, trie, bloomFilter, nodes);
	}

	private void checkTrie(File fromKraken, File reads, KrakenKMerFastqMerger krakenFilter, KMerTrie<String> trie,
			BloomFilter<byte[]> bloomFilter, Set<TaxIdNode> nodes) throws IOException {
		FileInputStream fStream2 = new FileInputStream(reads);
		GZIPInputStream gStream2 = new GZIPInputStream(fStream2, 4096);

		// Positive Test:
		krakenFilter.process(new BufferedInputStream(new FileInputStream(fromKraken)),
				new BufferedInputStream(gStream2),
				KrakenKMerFastqMerger.createFilterByTaxIdNodes(nodes, new FilterListener() {
					@Override
					public void newTaxidForRead(long readCount, String taxid, byte[] readDescriptor, byte[] read,
							byte[] readProbs) {
						assertEquals(taxid, trie.get(read, 0));
					}
				}));

		// Negative Test:
		byte[] read = new byte[trie.getLen()];
		byte[] cgat = { 'C', 'G', 'A', 'T' };
		Random random = new Random(42);

		for (int i = 1; i < 5 * 1000 * 1000; i++) {
			for (int j = 0; j < trie.getLen(); j++) {
				read[j] = cgat[random.nextInt(4)];
			}
			if (!bloomFilter.mightContain(read)) {
				assertNull(trie.get(read, 0));
			}
		}
	}
}
