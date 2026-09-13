/*
 * Genestrip
 */
package org.metagene.genestrip.refseq;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;
import org.metagene.genestrip.util.DigitTrie;

/**
 * Tests that a trie keyed by accessions accepts the bytes an accession is made of.
 * <p>
 * A {@link DigitTrie} maps a byte to {@code bite - '0'} over a range of ten, so a letter falls
 * outside it and the node lookup returns null - which a write then dereferences. Accession keys were
 * being stored in a plain {@code DigitTrie}, and the whole catalog pass died on the first finished
 * replicon it recognised.
 */
public class AccessionTrieTest {

	private static byte[] key(String s) {
		byte[] b = new byte[s.length()];
		for (int i = 0; i < s.length(); i++) {
			b[i] = (byte) s.charAt(i);
		}
		return b;
	}

	private static void put(AccessionTrie<byte[]> trie, String k, String v) {
		byte[] b = key(k);
		trie.set(b, 0, b.length, key(v));
	}

	private static byte[] get(AccessionTrie<byte[]> trie, String k) {
		byte[] b = key(k);
		return trie.get(b, 0, b.length);
	}

	/** The accessions a catalog actually carries must store and come back. */
	@Test
	public void testAccessionKeysRoundTrip() {
		AccessionTrie<byte[]> trie = new AccessionTrie<byte[]>();
		String[] keys = { "NC_001911", "NC_004843", "NZ_CABEIU01", "AC_000091", "NW_020955238", "NG_012345" };
		for (String k : keys) {
			put(trie, k, k);
		}
		for (String k : keys) {
			assertArrayEquals("key " + k + " must come back", key(k), get(trie, k));
		}
		assertNull("an absent key must be null", get(trie, "NC_999999"));
	}

	/** Keys sharing a prefix must stay distinct, which is the point of a trie. */
	@Test
	public void testSharedPrefixesStayDistinct() {
		AccessionTrie<byte[]> trie = new AccessionTrie<byte[]>();
		put(trie, "NC_0019", "short");
		put(trie, "NC_001911", "long");
		assertArrayEquals(key("short"), get(trie, "NC_0019"));
		assertArrayEquals(key("long"), get(trie, "NC_001911"));
	}

	/** No byte may be rejected: anything outside the alphabet shares a catch-all slot. */
	@Test
	public void testAnyByteIsAccepted() {
		AccessionTrie<byte[]> trie = new AccessionTrie<byte[]>();
		put(trie, "a.b|c-1", "odd");
		assertArrayEquals(key("odd"), get(trie, "a.b|c-1"));
	}

	/**
	 * The bug itself: a plain digit trie cannot take a letter, and fails on the write rather than
	 * refusing it. This is what the accession trie exists to avoid, so it is worth stating.
	 */
	@Test
	public void testAPlainDigitTrieCannotTakeAnAccession() {
		DigitTrie<byte[]> digits = new DigitTrie<byte[]>();
		byte[] b = key("NC_001911");
		try {
			digits.set(b, 0, b.length, key("x"));
			fail("a digit trie must not silently accept a letter key");
		} catch (NullPointerException expected) {
			assertTrue(true);
		}
		assertNull("and it stores nothing", digits.get(b, 0, b.length));
	}
}
