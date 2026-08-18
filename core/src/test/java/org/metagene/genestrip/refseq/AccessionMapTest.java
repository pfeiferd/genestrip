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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.function.Supplier;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameter;
import org.junit.runners.Parameterized.Parameters;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * Tests every {@link AccessionMap} implementation against the same expectations, so that the two
 * remain interchangeable: what a caller may rely on is the interface, not the layout behind it.
 */
@RunWith(Parameterized.class)
public class AccessionMapTest {
	/**
	 * The implementations under test.
	 *
	 * @return the name and the constructor of each implementation
	 */
	@Parameters(name = "{0}")
	@SuppressWarnings("deprecation")
	public static Collection<Object[]> implementations() {
		return Arrays.asList(new Object[][] {
				{ AccessionMapImpl.class.getSimpleName(), (Supplier<AccessionMap>) AccessionMapImpl::new },
				{ AccessionMapTrieImpl.class.getSimpleName(), (Supplier<AccessionMap>) AccessionMapTrieImpl::new } });
	}

	/** The name of the implementation under test, which names the test run. */
	@Parameter(0)
	public String name;

	/** Creates an empty map of the implementation under test. */
	@Parameter(1)
	public Supplier<AccessionMap> factory;

	private static TaxIdNode[] nodes;

	/**
	 * Returns three tax nodes to use as map values, from a tax tree written for the test.
	 *
	 * @return the nodes
	 * @throws IOException if the temporary tax tree cannot be written
	 */
	private static synchronized TaxIdNode[] nodes() throws IOException {
		if (nodes == null) {
			File dir = Files.createTempDirectory("accessionmaptest").toFile();
			dir.deleteOnExit();
			StringBuilder nodesDmp = new StringBuilder();
			StringBuilder namesDmp = new StringBuilder();
			for (int[] edge : new int[][] { { 1, 1 }, { 2, 1 }, { 3, 1 } }) {
				nodesDmp.append(edge[0]).append("\t|\t").append(edge[1]).append("\t|\tno rank\t|\t\t|\n");
				namesDmp.append(edge[0]).append("\t|\t").append(edge[0]).append("\t|\t\t|\tscientific name\t|\n");
			}
			Files.write(new File(dir, TaxTree.NODES_DMP).toPath(), nodesDmp.toString().getBytes(StandardCharsets.UTF_8));
			Files.write(new File(dir, TaxTree.NAMES_DMP).toPath(), namesDmp.toString().getBytes(StandardCharsets.UTF_8));
			TaxTree tree = new TaxTree(dir, false);
			nodes = new TaxIdNode[] { tree.getNodeByTaxId("1"), tree.getNodeByTaxId("2"), tree.getNodeByTaxId("3") };
		}
		return nodes;
	}

	private static byte[] bytes(String key) {
		return key.getBytes(StandardCharsets.UTF_8);
	}

	/**
	 * Fills a map with the given keys in the given order and returns what the map is then expected
	 * to answer. Which node a key gets follows from the key itself rather than from its position, so
	 * that filling the same keys in another order must yield the very same map.
	 */
	private Map<String, TaxIdNode> fill(AccessionMap map, List<String> keys) throws IOException {
		TaxIdNode[] values = nodes();
		Map<String, TaxIdNode> expected = new HashMap<>();
		for (String each : keys) {
			byte[] key = bytes(each);
			TaxIdNode node = values[Math.floorMod(each.hashCode(), values.length)];
			map.put(key, 0, key.length, node);
			expected.put(each, node);
		}
		map.optimize();
		return expected;
	}

	/** Generates distinct keys shaped like RefSeq accessions. */
	private static List<String> accessionLikeKeys(int count, long seed) {
		Random random = new Random(seed);
		String[] prefixes = { "NC_", "NZ_", "AC_", "NG_", "XM_", "NR_" };
		Set<String> keys = new LinkedHashSet<>();
		while (keys.size() < count) {
			StringBuilder key = new StringBuilder(prefixes[random.nextInt(prefixes.length)]);
			for (int i = 0; i < 6 + random.nextInt(6); i++) {
				key.append((char) ('0' + random.nextInt(10)));
			}
			if (random.nextInt(4) == 0) {
				key.append('.').append(random.nextInt(9) + 1);
			}
			keys.add(key.toString());
		}
		return new ArrayList<>(keys);
	}

	/**
	 * Every stored key must come back with the value it was stored with, whether it is looked up on
	 * its own or as a range inside a larger buffer.
	 *
	 * @throws IOException if the tax tree cannot be built
	 */
	@Test
	public void testStoredKeysAreFound() throws IOException {
		AccessionMap map = factory.get();
		Map<String, TaxIdNode> expected = fill(map, accessionLikeKeys(20000, 5));
		for (Map.Entry<String, TaxIdNode> entry : expected.entrySet()) {
			byte[] key = bytes(entry.getKey());
			assertSame(entry.getKey(), entry.getValue(), map.get(key, 0, key.length, false));
			byte[] embedded = bytes("prefix>" + entry.getKey() + "<suffix");
			assertSame(entry.getKey(), entry.getValue(), map.get(embedded, 7, 7 + entry.getKey().length(), false));
		}
	}

	/**
	 * Keys that were never stored must not be found, including ones that are a prefix of, an
	 * extension of, or a near miss for a stored key.
	 *
	 * @throws IOException if the tax tree cannot be built
	 */
	@Test
	public void testAbsentKeysAreNotFound() throws IOException {
		AccessionMap map = factory.get();
		List<String> keys = accessionLikeKeys(5000, 7);
		Map<String, TaxIdNode> expected = fill(map, keys);
		for (String stored : keys.subList(0, 200)) {
			String[] candidates = { stored + "0", stored.substring(0, stored.length() - 1),
					"Q" + stored.substring(1) };
			for (String candidate : candidates) {
				if (!candidate.isEmpty() && !expected.containsKey(candidate)) {
					byte[] key = bytes(candidate);
					assertNull(candidate, map.get(key, 0, key.length, false));
				}
			}
		}
	}

	/**
	 * The order in which entries are put must not change what the map answers.
	 *
	 * @throws IOException if the tax tree cannot be built
	 */
	@Test
	public void testPutOrderDoesNotMatter() throws IOException {
		List<String> keys = accessionLikeKeys(5000, 11);
		AccessionMap first = factory.get();
		Map<String, TaxIdNode> expected = fill(first, keys);
		List<String> shuffled = new ArrayList<>(keys);
		Collections.shuffle(shuffled, new Random(13));
		AccessionMap second = factory.get();
		fill(second, shuffled);
		for (String stored : keys) {
			byte[] key = bytes(stored);
			assertSame(stored, expected.get(stored), second.get(key, 0, key.length, false));
		}
	}

	/**
	 * Keys shorter than any prefix an implementation may group by, or holding characters other than
	 * the ones accessions are made of, must be stored and found like any other.
	 *
	 * @throws IOException if the tax tree cannot be built
	 */
	@Test
	public void testUnusualKeysAreStored() throws IOException {
		AccessionMap map = factory.get();
		List<String> keys = new ArrayList<>(Arrays.asList("A", "AB", "ABC", "ABCD", "ABCDE", "ABCDEF", "lower_case",
				"has space", "tab\tinside", "ümlaut", "0", "..", "__"));
		keys.addAll(accessionLikeKeys(500, 17));
		Map<String, TaxIdNode> expected = fill(map, keys);
		for (Map.Entry<String, TaxIdNode> entry : expected.entrySet()) {
			byte[] key = bytes(entry.getKey());
			assertSame(entry.getKey(), entry.getValue(), map.get(key, 0, key.length, false));
		}
	}

	/**
	 * Keys too long for whatever slot an implementation keeps them in must still be ordered and
	 * found correctly, including ones that agree far enough that only the whole key tells them apart.
	 *
	 * @throws IOException if the tax tree cannot be built
	 */
	@Test
	public void testLongKeysAreStored() throws IOException {
		Random random = new Random(19);
		Set<String> keys = new LinkedHashSet<>();
		for (int length : new int[] { 16, 17, 24, 25, 40, 120, 255 }) {
			for (int i = 0; i < 200; i++) {
				StringBuilder key = new StringBuilder("NC_");
				while (key.length() < length) {
					key.append((char) ('A' + random.nextInt(26)));
				}
				keys.add(key.substring(0, length));
			}
		}
		StringBuilder shared = new StringBuilder("NC_");
		while (shared.length() < 60) {
			shared.append('Q');
		}
		for (int i = 1000; i < 1200; i++) {
			keys.add(shared.toString() + i);
		}
		AccessionMap map = factory.get();
		Map<String, TaxIdNode> expected = fill(map, new ArrayList<>(keys));
		for (Map.Entry<String, TaxIdNode> entry : expected.entrySet()) {
			byte[] key = bytes(entry.getKey());
			assertSame(entry.getKey(), entry.getValue(), map.get(key, 0, key.length, false));
		}
	}

	/**
	 * A key of the greatest supported length must be stored, and a longer one must be refused rather
	 * than kept in a way that could not be found again.
	 *
	 * @throws IOException if the tax tree cannot be built
	 */
	@Test
	public void testOverlongKeyIsRefused() throws IOException {
		StringBuilder longest = new StringBuilder("NC_");
		while (longest.length() < 255) {
			longest.append('Y');
		}
		AccessionMap map = factory.get();
		byte[] key = bytes(longest.toString());
		map.put(key, 0, key.length, nodes()[0]);
		map.optimize();
		assertSame(nodes()[0], map.get(key, 0, key.length, false));

		StringBuilder tooLong = new StringBuilder("NC_");
		while (tooLong.length() < 400) {
			tooLong.append('Z');
		}
		byte[] refused = bytes(tooLong.toString());
		try {
			factory.get().put(refused, 0, refused.length, nodes()[0]);
			fail("a key of " + refused.length + " bytes should have been refused");
		} catch (IllegalArgumentException expected) {
			assertTrue(expected.getMessage(), expected.getMessage().contains("length"));
		}
	}

	/**
	 * Looking up before the map is optimized must fail rather than quietly answer from entries that
	 * are not in order yet.
	 *
	 * @throws IOException if the tax tree cannot be built
	 */
	@Test
	public void testLookupBeforeOptimizeFails() throws IOException {
		AccessionMap map = factory.get();
		byte[] key = bytes("NC_012345");
		map.put(key, 0, key.length, nodes()[0]);
		try {
			map.get(key, 0, key.length, false);
			fail("a lookup before optimize() should have failed");
		} catch (IllegalStateException expected) {
			assertNotNull(expected.getMessage());
		}
	}

	/**
	 * The entries of a node must be counted across the whole map.
	 *
	 * @throws IOException if the tax tree cannot be built
	 */
	@Test
	public void testEntriesPerNodeAreCounted() throws IOException {
		AccessionMap map = factory.get();
		List<String> keys = accessionLikeKeys(3000, 23);
		Map<String, TaxIdNode> expected = fill(map, keys);
		int total = 0;
		for (TaxIdNode node : nodes()) {
			int counted = 0;
			for (TaxIdNode assigned : expected.values()) {
				if (assigned == node) {
					counted++;
				}
			}
			assertTrue("every node should hold entries", counted > 0);
			assertEquals(node.getTaxId(), counted, map.getEntriesForNode(node));
			total += map.getEntriesForNode(node);
		}
		assertEquals(keys.size(), total);
	}

	/**
	 * With {@code assemblyAccessionsOnly} set, only accessions whose prefix marks a complete genomic,
	 * an RNA or a messenger-RNA sequence may be answered.
	 *
	 * @throws IOException if the tax tree cannot be built
	 */
	@Test
	public void testCompleteGenomesOnlyFiltersByPrefix() throws IOException {
		AccessionMap map = factory.get();
		fill(map, Arrays.asList("NC_012345", "NZ_012345", "AC_012345", "NR_012345", "NM_012345", "NG_012345",
				"XP_012345"));
		for (String accepted : new String[] { "NC_012345", "NZ_012345", "AC_012345", "NR_012345", "NM_012345" }) {
			byte[] key = bytes(accepted);
			assertNotNull(accepted, map.get(key, 0, key.length, true));
		}
		for (String filtered : new String[] { "NG_012345", "XP_012345" }) {
			byte[] key = bytes(filtered);
			assertNull(filtered, map.get(key, 0, key.length, true));
			assertNotNull(filtered, map.get(key, 0, key.length, false));
		}
	}
}
