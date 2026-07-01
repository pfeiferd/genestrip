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

/**
 * A {@link DigitTrie} that maps digit keys to {@link StringLong} counters, providing increment and
 * add operations that create the counter on first use.
 */
public class StringLongDigitTrie extends DigitTrie<StringLong> {
	/**
	 * Creates an empty trie.
	 */
	public StringLongDigitTrie() {
	}

	/**
	 * Increments the counter for the key given by {@code seq[start, end)} by one, creating it if
	 * necessary.
	 *
	 * @param seq the byte array containing the key.
	 * @param start the start index of the key (inclusive).
	 * @param end the end index of the key (exclusive).
	 * @return the (possibly newly created) counter after incrementing.
	 */
	public StringLong inc(byte[] seq, int start, int end) {
		return add(seq, start, end, 1);
	}

	/**
	 * Increments the counter for the given key by one, creating it if necessary.
	 *
	 * @param key the digit key.
	 * @return the (possibly newly created) counter after incrementing.
	 */
	public StringLong inc(String key) {
		return add(key, 1);
	}
	
	/**
	 * Adds {@code add} to the counter for the key given by {@code seq[start, end)}, creating it if
	 * necessary.
	 *
	 * @param seq the byte array containing the key.
	 * @param start the start index of the key (inclusive).
	 * @param end the end index of the key (exclusive).
	 * @param add the amount to add to the counter.
	 * @return the (possibly newly created) counter after adding.
	 */
	public StringLong add(byte[] seq, int start, int end, long add) {
		StringLong stringLong = get(seq, start, end, this);
		stringLong.longValue += add;

		return stringLong;
	}

	/**
	 * Adds {@code add} to the counter for the given key, creating it if necessary.
	 *
	 * @param key the digit key.
	 * @param add the amount to add to the counter.
	 * @return the (possibly newly created) counter after adding.
	 */
	public StringLong add(String key, long add) {
		StringLong stringLong = get(key, this);
		stringLong.longValue += add;

		return stringLong;
	}
	
	@Override
	protected StringLong createInGet(byte[] seq, int start, int end, Object createContext) {
		return new StringLong(new String(seq, start, end - start));
	}

	@Override
	protected StringLong createInGet(String digits, Object createContext) {
		return new StringLong(digits);
	}

	/**
	 * A mutable pairing of a string key and an associated long value (counter).
	 */
	public static class StringLong implements Serializable, Comparable<StringLong> {
		private static final long serialVersionUID = 1L;
		/** The string key of this pairing. */
		protected String stringValue;
		/** The long value (counter) associated with the key. */
		protected long longValue;

		/**
		 * Creates a pairing for the given string key with a zero value.
		 *
		 * @param stringValue the string key.
		 */
		public StringLong(String stringValue) {
			this.stringValue = stringValue;
		}

		/**
		 * Returns the long value (counter) of this pairing.
		 *
		 * @return the long value.
		 */
		public long getLongValue() {
			return longValue;
		}

		/**
		 * Returns the string key of this pairing.
		 *
		 * @return the string key.
		 */
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
