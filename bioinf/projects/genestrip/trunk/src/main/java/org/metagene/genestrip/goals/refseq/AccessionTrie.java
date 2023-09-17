package org.metagene.genestrip.goals.refseq;

import java.io.Serializable;

import org.metagene.genestrip.util.DigitTrie;

public class AccessionTrie<V extends Serializable> extends DigitTrie<V> {
	private static final long serialVersionUID = 1L;
	
	@Override
	protected int mapToIndex(byte bite, int pos) {
		if (pos <= 1) {
			return bite - 'A';
		}
		if (pos == 2) {
			return bite - '_';
		}
		return bite == '.' ? 10 : bite - '0';
	}

	@Override
	protected int range(int pos) {
		if (pos <= 1) {
			return 26;
		}
		if (pos == 2) {
			return 1;
		}
		return 11;
	}
}
