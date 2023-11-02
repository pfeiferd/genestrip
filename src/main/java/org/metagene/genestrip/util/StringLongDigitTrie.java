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
package org.metagene.genestrip.util;

import java.io.Serializable;

import org.metagene.genestrip.util.StringLongDigitTrie.StringLong;

public class StringLongDigitTrie extends DigitTrie<StringLong> {
	public StringLong inc(byte[] seq, int start, int end) {
		return add(seq, start, end, 1);
	}

	public StringLong inc(String key) {
		return add(key, 1);
	}
	
	public StringLong add(byte[] seq, int start, int end, long add) {
		StringLong stringLong = get(seq, start, end, this);
		stringLong.longValue += add;

		return stringLong;
	}

	public StringLong add(String key, long add) {
		StringLong stringLong = get(key, this);
		stringLong.longValue += add;

		return stringLong;
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

	public static class StringLong implements Serializable, Comparable<StringLong> {
		private static final long serialVersionUID = 1L;
		private String stringValue;
		private long longValue;

		public long getLongValue() {
			return longValue;
		}

		public String getStringValue() {
			return stringValue;
		}
		
		@Override
		public int compareTo(StringLong o) {
			if (stringValue == null) {
				return -1;
			}
			if (o.stringValue == null) {
				return 1;
			}
			return stringValue.compareTo(o.stringValue);
		}
		
		@Override
		public String toString() {
			return "SL:(" + stringValue + ", " + longValue + ")";
		}
	}
}
