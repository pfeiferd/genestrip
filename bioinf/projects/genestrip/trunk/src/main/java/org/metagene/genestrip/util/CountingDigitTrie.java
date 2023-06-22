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

import java.util.Map;

public class CountingDigitTrie {
	private CountingDigitTrie[] children;
	private long value;
	private String digits;

	public CountingDigitTrie() {
	}

	public String inc(String digits) {
		return add(digits, 1);
	}

	public String add(String digits, int add) {
		int index;
		int end = digits.length();
		CountingDigitTrie node = this, child;
		for (int i = 0; i < end; i++, node = child) {
			index = digits.charAt(i) - '0';
			if (index > 9 || index < 0) {
				return null;
			}
			if (node.children == null) {
				node.children = new CountingDigitTrie[10];
			}
			child = node.children[index];
			if (child == null) {
				child = new CountingDigitTrie();
				node.children[index] = child;
			}
		}
		node.value += add;
		if (node.digits == null) {
			node.digits = digits;
		}

		return node.digits;
	}

	public String inc(byte[] seq, int start, int end) {
		return add(seq, start, end, 1);
	}

	public String add(byte[] seq, int start, int end, int add) {
		int index;
		CountingDigitTrie node = this, child;
		for (int i = start; i < end; i++, node = child) {
			index = seq[i] - '0';
			if (index > 9 || index < 0) {
				return null;
			}
			if (node.children == null) {
				node.children = new CountingDigitTrie[10];
			}
			child = node.children[index];
			if (child == null) {
				child = new CountingDigitTrie();
				node.children[index] = child;
			}
		}
		node.value += add;
		if (node.digits == null) {
			node.digits = new String(seq, start, end - start);
		}
		return node.digits;
	}

	public String get(byte[] seq, int start, int end) {
		int index;
		CountingDigitTrie node = this, child;
		for (int i = start; i < end; i++, node = child) {
			index = seq[i] - '0';
			if (index > 9 || index < 0) {
				return null;
			}
			if (node.children == null) {
				node.children = new CountingDigitTrie[10];
			}
			child = node.children[index];
			if (child == null) {
				child = new CountingDigitTrie();
				node.children[index] = child;
			}
		}
		if (node.digits == null) {
			node.digits = new String(seq, start, end - start);
		}
		return node.digits;
	}

	public void collect(Map<String, Long> map) {
		if (digits != null) {
			map.put(digits, value);
		}
		if (children != null) {
			for (int k = 0; k < children.length; k++) {
				if (children[k] != null) {
					children[k].collect(map);
				}
			}
		}
	}
}
