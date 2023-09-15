package org.metagene.genestrip.goals.refseq;

import java.io.Serializable;

import org.metagene.genestrip.util.DigitTrie;

public class AccessionTrie<V extends Serializable> extends DigitTrie<V> {
	private static final long serialVersionUID = 1L;
	
	private static int[] JUMP_TABLE = new int[128];
	private static int RANGE;
	{
		RANGE = 0;
		for (int i = 0; i < JUMP_TABLE.length; i++) {
			if ((Character.isAlphabetic(i) && Character.isUpperCase(i)) || i == '.' || i == '_'
					|| Character.isDigit(i)) {
				JUMP_TABLE[i] = RANGE;
				RANGE++;
			}
		}
	}

	@Override
	protected int mapToIndex(byte bite) {
		return JUMP_TABLE[bite];
	}

	@Override
	protected int range() {
		return RANGE;
	}
}
