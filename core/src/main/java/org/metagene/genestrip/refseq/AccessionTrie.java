/*
 * Genestrip
 */
package org.metagene.genestrip.refseq;

import org.metagene.genestrip.util.DigitTrie;

/**
 * A trie keyed by RefSeq accessions rather than by digits alone.
 * <p>
 * {@link DigitTrie} maps a byte to {@code bite - '0'} over a range of ten, so a letter lands outside
 * the range and the lookup yields nothing - for a write, a null node that the caller then
 * dereferences. An accession is letters, digits and an underscore, so it needs an alphabet of its
 * own, which is what this adds. Anything outside that alphabet shares one catch-all slot, so no byte
 * is ever rejected and no key can fail to be stored.
 *
 * @param <V> the type of the values stored.
 */
public class AccessionTrie<V> extends DigitTrie<V> {

	/** Creates an empty trie. */
	public AccessionTrie() {
	}

	/** Digits, then the upper-case letters, then the underscore, then anything else. */
	@Override
	protected int mapToIndex(byte bite, int pos) {
		if (bite >= '0' && bite <= '9') {
			return bite - '0';
		}
		if (bite >= 'A' && bite <= 'Z') {
			return bite - 'A' + 10;
		}
		return bite == '_' ? 36 : 37;
	}

	@Override
	protected int range(int pos) {
		return 38;
	}
}
