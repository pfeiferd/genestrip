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
import java.util.Collection;

public class DigitTrie<V extends Serializable> implements Serializable {
	private static final long serialVersionUID = 1L;
	
	private DigitTrie<V>[] children;
	private V value;

	public DigitTrie() {
	}

	public V set(String digits, V value) {
		DigitTrie<V> node = getNode(digits, true);
		V oldValue = node.value;
		node.value = value;

		return oldValue;
	}

	public V set(byte[] seq, int start, int end, V value) {
		DigitTrie<V> node = getNode(seq, start, end, true);
		V oldValue = node.value;
		node.value = value;

		return oldValue;
	}

	public V get(byte[] seq, int start, int end) {
		return get(seq, start, end, false);
	}

	public V get(byte[] seq, int start, int end, boolean create) {
		DigitTrie<V> node = getNode(seq, start, end, create);
		if (node == null) {
			if (create) {
				throw new IllegalStateException("Cant get node via sequence " + new String(seq, start, end - start));				
			}
			return null;
		}
		if (node.value == null && create) {
			node.value = createInGet(seq, start, end);
		}

		return node.value;
	}

	@SuppressWarnings("unchecked")
	private DigitTrie<V> getNode(byte[] seq, int start, int end, boolean create) {
		int index;
		DigitTrie<V> node = this, child;
		for (int i = start; i < end; i++, node = child) {
			index = mapToIndex(seq[i]);
			if (index < 0 || index >= range()) {
				return null;
			}
			if (node.children == null) {
				if (create) {
					node.children = new DigitTrie[range()];
				} else {
					return null;
				}
			}
			child = node.children[index];
			if (child == null) {
				if (create) {
					child = new DigitTrie<V>();
					node.children[index] = child;
				} else {
					return null;
				}
			}
		}
		return node;
	}
	
	protected int mapToIndex(byte bite) {
		return bite - '0';
	}
	
	protected int range() {
		return 10;
	}

	protected V createInGet(byte[] seq, int start, int end) {
		return null;
	}

	public V get(String digits) {
		return get(digits, false);
	}

	public V get(String digits, boolean create) {
		DigitTrie<V> node = getNode(digits, create);
		if (node == null) {
			return null;
		}
		if (node.value == null && create) {
			node.value = createInGet(digits);
		}

		return node.value;
	}

	@SuppressWarnings("unchecked")
	private DigitTrie<V> getNode(String digits, boolean create) {
		int index;
		int end = digits.length();
		DigitTrie<V> node = this, child;
		for (int i = 0; i < end; i++, node = child) {
			index = mapToIndex((byte) digits.charAt(i));
			if (index < 0 || index >= range()) {
				return null;
			}
			if (node.children == null) {
				if (create) {
					node.children = new DigitTrie[range()];
				} else {
					return null;
				}
			}
			child = node.children[index];
			if (child == null) {
				if (create) {
					child = new DigitTrie<V>();
					node.children[index] = child;
				} else {
					return null;
				}
			}
		}
		return node;
	}

	protected V createInGet(String digits) {
		return null;
	}

	public void collect(Collection<V> collection) {
		if (value != null) {
			collection.add(value);
		}
		if (children != null) {
			for (int k = 0; k < children.length; k++) {
				if (children[k] != null) {
					children[k].collect(collection);
				}
			}
		}
	}
}
