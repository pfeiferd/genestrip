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
package org.metagene.genestrip.trie;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.metagene.genestrip.util.CGAT;

public class KMerTrie<V extends Serializable> implements Serializable {
	private static final long serialVersionUID = 1L;

	private final int factor;
	private final int len;
	private final Object[] root;
	private final int[] jumpTable;
	private final boolean allowDoubleEntry;
	private boolean compressed;
	private long entries;

	public KMerTrie(int len) {
		this(2, len, false);
	}

	public KMerTrie(int factor, int len, boolean allowDoubleEntry) {
		if (factor > 3 || factor < 1) {
			throw new IllegalArgumentException("factor must be >= 1 and <= 3");
		}
		jumpTable = new int[Byte.MAX_VALUE];
		for (int i = 0; i < jumpTable.length; i++) {
			jumpTable[i] = -1;
		}
		jumpTable['C'] = 0;
		jumpTable['G'] = 1;
		jumpTable['A'] = 2;
		jumpTable['T'] = 3;

		this.factor = factor;
		this.len = len;
		this.allowDoubleEntry = allowDoubleEntry;

		int size = 1;
		for (int i = 0; i < factor; i++) {
			size *= 4;
		}
		root = new Object[size];
	}

	public int getLen() {
		return len;
	}

	public long getEntries() {
		return entries;
	}

	public void put(byte[] nseq, int start, V value) {
		if (compressed) {
			throw new IllegalStateException("Cant insert in compressed trie");
		}

		Object[] node = root;

		int mult;
		int pos = 0;
		int j = 0;

		for (int i = 0; i < len; i += factor) {
			pos = 0;
			mult = 1;
			for (j = 0; j < factor && i + j < len; j++) {
				byte c = CGAT.cgatToUpperCase(nseq[start + i + j]);
				if (c < 0 || jumpTable[c] == -1) {
					throw new IllegalArgumentException("Not a CGAT sequence");
				}
				pos += jumpTable[c] * mult;
				mult *= 4;
			}
			Object[] next = (Object[]) node[pos];
			if (next == null) {
				next = new Object[root.length];
				node[pos] = next;
			}
			node = next;
		}
		if (!allowDoubleEntry && node[pos] != null) {
			throw new IllegalStateException("double entry " + node[pos] + " " + value);
		}
		if (node[pos] == null) {
			entries++;
		}
		node[pos] = value;
	}

	@SuppressWarnings("unchecked")
	public V get(byte[] nseq, int start) {
		Object node = root;

		int mult;
		int pos = 0;
		int j = 0;
		int snPosIndex = 0;

		for (int i = 0; i < len && node != null; i += factor) {
			pos = 0;
			mult = 1;
			for (j = 0; j < factor && i + j < len; j++) {
				byte c = CGAT.cgatToUpperCase(nseq[start + i + j]);
				if (c < 0 || jumpTable[c] == -1) {
					return null;
				}
				pos += jumpTable[c] * mult;
				mult *= 4;
			}
			if (node instanceof Object[]) {
				node = ((Object[]) node)[pos];
			} else {
				SmallNode sn = (SmallNode) node;
				node = null;
				if (sn.getPosVal(snPosIndex) == pos) {
					snPosIndex = sn.nextPos(snPosIndex);
					if (snPosIndex == 0) {
						node = sn.next;
					} else {
						node = sn;
					}
				}
			}
		}
		if (node == null) {
			return null;
		} else if (node instanceof Object[]) {
			return (V) ((Object[]) node)[pos];
		} else {
			return (V) ((SmallNode) node).next;
		}
	}

	public void compress() {
		compressHelp(root);
		compressed = true;
	}

	private void compressHelp(Object[] parent) {
		for (int j = 0; j < parent.length; j++) {
			if (parent[j] instanceof Object[]) {
				Object[] n = (Object[]) parent[j];
				compressHelp(n);

				int notNulls = 0;
				int notNullPos = 0;
				for (int i = 0; i < n.length; i++) {
					if (n[i] != null) {
						notNulls++;
						if (notNulls == 1) {
							notNullPos = i;
						} else if (notNulls > 1) {
							break;
						}
					}
				}
				if (notNulls == 0) {
					throw new IllegalStateException("There should be at least one successor");
				} else if (notNulls == 1) {
					SmallNode sn;
					if (n[notNullPos] instanceof SmallNode) {
						sn = (SmallNode) n[notNullPos];
						sn = sn.add(notNullPos);
					} else {
						sn = new SmallNode(notNullPos, n[notNullPos]);
					}
					parent[j] = sn;
				}
			}
		}
	}

	public static class SmallNode implements Serializable {
		private static final long serialVersionUID = 1L;

		private byte pos0;
		private byte pos1;
		private byte pos2;
		private byte pos3;
		private byte pos4;
		private byte pos5;
		private byte pos6;
		private byte pos7;
		private Object next;

		public SmallNode(int pos, Object next) {
			pos0 = (byte) pos;
			pos1 = pos2 = pos3 = pos4 = pos5 = pos6 = pos7 = -1;
			this.next = next;
		}

		public SmallNode add(int pos) {
			if (pos1 == -1) {
				pos1 = (byte) pos0;
				pos0 = (byte) pos;
				return this;
			}
			if (pos2 == -1) {
				pos2 = (byte) pos1;
				pos1 = (byte) pos0;
				pos0 = (byte) pos;
				return this;
			}
			if (pos3 == -1) {
				pos3 = (byte) pos2;
				pos2 = (byte) pos1;
				pos1 = (byte) pos0;
				pos0 = (byte) pos;
				return this;
			}
			if (pos4 == -1) {
				pos4 = (byte) pos3;
				pos3 = (byte) pos2;
				pos2 = (byte) pos1;
				pos1 = (byte) pos0;
				pos0 = (byte) pos;
				return this;
			}
			if (pos5 == -1) {
				pos5 = (byte) pos4;
				pos4 = (byte) pos3;
				pos3 = (byte) pos2;
				pos2 = (byte) pos1;
				pos1 = (byte) pos0;
				pos0 = (byte) pos;
				return this;
			}
			if (pos6 == -1) {
				pos6 = (byte) pos5;
				pos5 = (byte) pos4;
				pos4 = (byte) pos3;
				pos3 = (byte) pos2;
				pos2 = (byte) pos1;
				pos1 = (byte) pos0;
				pos0 = (byte) pos;
				return this;
			}
			if (pos7 == -1) {
				pos7 = (byte) pos6;
				pos6 = (byte) pos5;
				pos5 = (byte) pos4;
				pos4 = (byte) pos3;
				pos3 = (byte) pos2;
				pos2 = (byte) pos1;
				pos1 = (byte) pos0;
				pos0 = (byte) pos;
				return this;
			}
			return new SmallNode(pos, this);
		}

		public byte getPosVal(int posIndex) {
			if (posIndex == 0) {
				return pos0;
			}
			if (posIndex == 1) {
				return pos1;
			}
			if (posIndex == 2) {
				return pos2;
			}
			if (posIndex == 3) {
				return pos3;
			}
			if (posIndex == 4) {
				return pos4;
			}
			if (posIndex == 5) {
				return pos5;
			}
			if (posIndex == 6) {
				return pos6;
			}
			if (posIndex == 7) {
				return pos7;
			}
			throw new IllegalArgumentException("Inavlid pos index");
		}

		public int nextPos(int posIndex) {
			int next = (posIndex + 1) % 8;
			return getPosVal(next) == -1 ? 0 : next;
		}
	}

	public void writeAsJSON(PrintStream out) {
		writeAsJSONHelp(out, root);
	}

	protected void writeAsJSONHelp(PrintStream out, Object node) {
		if (node instanceof Object[]) {
			Object[] array = (Object[]) node;
			out.print('[');
			for (int i = 0; i < array.length; i++) {
				if (i != 0) {
					out.print(',');
				}
				writeAsJSONHelp(out, array[i]);
			}
			out.print(']');
		} else if (node instanceof SmallNode) {
			out.print('[');
			out.print('\"');
			SmallNode sn = (SmallNode) node;
			for (int pos = 0; pos != -1; pos = sn.nextPos(pos)) {
				out.print((char) sn.getPosVal(pos));
			}
			out.print('\"');
			out.print(',');
			writeAsJSONHelp(out, sn.next);
			out.print(']');
		} else if (node instanceof String) {
			out.print('\"');
			out.print(node);
			out.print('\"');
		} else {
			valueToJSON(out, node);
		}
	}
	
	protected void valueToJSON(PrintStream out, Object value) {
		out.print(value);
	}

	public void save(File trieFile) throws IOException {
		ObjectOutputStream oOut = new ObjectOutputStream(
				new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(trieFile))));
		oOut.writeObject(this);
		oOut.close();
	}

	public static KMerTrie<?> load(File filterFile) throws IOException, ClassNotFoundException {
		ObjectInputStream oOut = new ObjectInputStream(
				new BufferedInputStream(new GZIPInputStream(new FileInputStream(filterFile))));
		KMerTrie<?> res = (KMerTrie<?>) oOut.readObject();
		oOut.close();
		return res;
	}
}