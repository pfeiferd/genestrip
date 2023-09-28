package org.metagene.genestrip.kraken;

import java.io.Serializable;

import org.metagene.genestrip.kraken.StringLongDigitTrie.StringLong;
import org.metagene.genestrip.util.DigitTrie;

public class StringLongDigitTrie extends DigitTrie<StringLong> {
	public String inc(byte[] seq, int start, int end) {
		return add(seq, start, end, 1);
	}

	public String add(byte[] seq, int start, int end, long add) {
		StringLong stringLong = get(seq, start, end, true);
		stringLong.longValue += add;

		return stringLong.stringValue;
	}

	@Override
	protected StringLong createInGet(byte[] seq, int start, int end, Object createContext) {
		StringLong stringLong = new StringLong();
		stringLong.stringValue = new String(seq, start, end - start);

		return stringLong;
	}

	@Override
	protected StringLong createInGet(String digits, Object createContext) {
		StringLong stringLong = new StringLong();
		stringLong.stringValue = digits;

		return stringLong;
	}

	public static class StringLong implements Serializable {
		private static final long serialVersionUID = 1L;
		private String stringValue;
		private long longValue;

		public long getLongValue() {
			return longValue;
		}

		public String getStringValue() {
			return stringValue;
		}
	}
}
