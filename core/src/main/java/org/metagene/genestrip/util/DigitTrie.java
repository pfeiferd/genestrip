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

/**
 * A trie keyed on sequences of digit characters (0-9 by default) that maps each key to a value.
 * Subclasses may override {@link #mapToIndex(byte, int)}, {@link #range(int)} and the
 * {@code createInGet} hooks to support other alphabets and lazy value creation.
 *
 * @param <V> the type of the stored values.
 */
public class DigitTrie<V> implements Serializable {
	private static final long serialVersionUID = 1L;

	/** The child nodes indexed by mapped digit position, or null if this node has no children yet. */
	private DigitTrie<V>[] children;
	/** The value stored at this node, or null if no value is associated with the corresponding key. */
	private V value;

	/** Creates an empty trie with no stored values. */
	public DigitTrie() {
	}

	/**
	 * Associates the given value with the given digit key, returning the previously stored value or
	 * null.
	 *
	 * @param digits the digit key to associate the value with.
	 * @param value the value to store.
	 * @return the value previously stored for the key, or null if none.
	 */
	public V set(String digits, V value) {
		DigitTrie<V> node = getNode(digits, true);
		V oldValue = node.value;
		node.value = value;

		return oldValue;
	}

	/**
	 * Associates the given value with the key given by {@code seq[start, end)}, returning the
	 * previously stored value or null.
	 *
	 * @param seq the byte array containing the key.
	 * @param start the start index of the key (inclusive).
	 * @param end the end index of the key (exclusive).
	 * @param value the value to store.
	 * @return the value previously stored for the key, or null if none.
	 */
	public V set(byte[] seq, int start, int end, V value) {
		DigitTrie<V> node = getNode(seq, start, end, true);
		V oldValue = node.value;
		node.value = value;

		return oldValue;
	}

	/**
	 * Returns the value for the key given by {@code seq[start, end)}, or null if absent.
	 *
	 * @param seq the byte array containing the key.
	 * @param start the start index of the key (inclusive).
	 * @param end the end index of the key (exclusive).
	 * @return the value stored for the key, or null if absent.
	 */
	public V get(byte[] seq, int start, int end) {
		return get(seq, start, end, null);
	}

	/**
	 * Returns the value for the key given by {@code seq[start, end)}. If {@code createContext} is
	 * non-null and no value is stored yet, one is lazily created via
	 * {@link #createInGet(byte[], int, int, Object)} and stored.
	 *
	 * @param seq the byte array containing the key.
	 * @param start the start index of the key (inclusive).
	 * @param end the end index of the key (exclusive).
	 * @param createContext the context passed to lazy value creation, or null to disable creation.
	 * @return the stored or lazily created value, or null if absent and no creation was requested.
	 */
	public V get(byte[] seq, int start, int end, Object createContext) {
		DigitTrie<V> node = getNode(seq, start, end, createContext != null);
		if (node == null) {
			if (createContext != null) {
				throw new IllegalStateException("Cant get node via sequence " + new String(seq, start, end - start));				
			}
			return null;
		}
		if (node.value == null && createContext != null) {
			synchronized (node) {
				if (node.value == null) {
					node.value = createInGet(seq, start, end, createContext);
				}
			}
		}

		return node.value;
	}

	/**
	 * Returns the value for the given digit key. If {@code createContext} is non-null and no value is
	 * stored yet, one is lazily created via {@link #createInGet(String, Object)} and stored.
	 *
	 * @param digits the digit key to look up.
	 * @param createContext the context passed to lazy value creation, or null to disable creation.
	 * @return the stored or lazily created value, or null if absent and no creation was requested.
	 */
	public V get(String digits, Object createContext) {
		if (digits == null) {
			return null;
		}
		DigitTrie<V> node = getNode(digits, createContext != null);
		if (node == null) {
			return null;
		}
		if (node.value == null && createContext != null) {
			synchronized (node) {
				if (node.value == null) {
					node.value = createInGet(digits, createContext);
				}
			}
		}

		return node.value;
	}

	@SuppressWarnings("unchecked")
	private DigitTrie<V> getNode(byte[] seq, int start, int end, boolean create) {
		int index, pos;
		DigitTrie<V> node = this, child;
		for (int i = start; i < end; i++, node = child) {
			pos = i - start;
			index = mapToIndex(seq[i], pos);
			if (index < 0 || index >= range(pos)) {
				return null;
			}
			if (node.children == null) {
				if (create) {
					synchronized (node) {
						if (node.children == null) {
							node.children = new DigitTrie[range(pos)];
						}
					}
				} else {
					return null;
				}
			}
			child = node.children[index];
			if (child == null) {
				if (create) {
					synchronized (node) {
						child = node.children[index];
						if (child == null) {
							child = new DigitTrie<V>();
							node.children[index] = child;
						}
					}
				} else {
					return null;
				}
			}
		}
		return node;
	}

	@SuppressWarnings("unchecked")
	private DigitTrie<V> getNode(String digits, boolean create) {
		int index;
		int end = digits.length();
		DigitTrie<V> node = this, child;
		for (int i = 0; i < end; i++, node = child) {
			index = mapToIndex((byte) digits.charAt(i), i);
			if (index < 0 || index >= range(i)) {
				return null;
			}
			if (node.children == null) {
				if (create) {
					synchronized (node) {
						if (node.children == null) {
							node.children = new DigitTrie[range(i)];
						}
					}
				} else {
					return null;
				}
			}
			child = node.children[index];
			if (child == null) {
				if (create) {
					synchronized (node) {
						child = node.children[index];
						if (child == null) {
							child = new DigitTrie<V>();
							node.children[index] = child;
						}
					}
				} else {
					return null;
				}
			}
		}
		return node;
	}

	/**
	 * Maps a key byte at the given position to a child-array index; by default the digit value
	 * {@code bite - '0'}.
	 *
	 * @param bite the key byte to map.
	 * @param pos the position of the byte within the key.
	 * @return the child-array index for the byte.
	 */
	protected int mapToIndex(byte bite, int pos) {
		return bite - '0';
	}
	
	/**
	 * Returns the number of possible child indices at the given position; by default 10 (decimal
	 * digits).
	 *
	 * @param pos the position within the key.
	 * @return the number of possible child indices at that position.
	 */
	protected int range(int pos) {
		return 10;
	}

	/**
	 * Hook that creates the value to store when {@link #get(byte[], int, int, Object)} is called with
	 * a non-null create context and no value exists yet; returns null by default.
	 *
	 * @param seq the byte array containing the key.
	 * @param start the start index of the key (inclusive).
	 * @param end the end index of the key (exclusive).
	 * @param createContext the context passed in from the {@code get} call.
	 * @return the value to store, or null to store none.
	 */
	protected V createInGet(byte[] seq, int start, int end, Object createContext) {
		return null;
	}

	/**
	 * Returns the value for the given digit key, or null if absent.
	 *
	 * @param digits the digit key to look up.
	 * @return the value stored for the key, or null if absent.
	 */
	public V get(String digits) {
		return get(digits, null);
	}

	/**
	 * Hook that creates the value to store when {@link #get(String, Object)} is called with a
	 * non-null create context and no value exists yet; returns null by default.
	 *
	 * @param digits the digit key being looked up.
	 * @param createContext the context passed in from the {@code get} call.
	 * @return the value to store, or null to store none.
	 */
	protected V createInGet(String digits, Object createContext) {
		return null;
	}

	/**
	 * Adds all values stored in this trie (including descendants) to the given collection.
	 *
	 * @param collection the collection to add the stored values to.
	 */
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
