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

import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.io.Serializable;

import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATRingBuffer;
import org.metagene.genestrip.util.StreamProvider;

public class KMerTrie<V extends Serializable> implements Serializable {
	private static final InternalNullMarker NULL = new InternalNullMarker();

	private static final long serialVersionUID = 1L;

	private static final int[] CGAT_JUMP_TABLE;
	private static final int[] CGAT_REVERSE_JUMP_TABLE;
	private static final byte[] DECODE_TABLE, REVERSE_DECODE_TABLE;

	static {
		CGAT_JUMP_TABLE = new int[Byte.MAX_VALUE];
		CGAT_REVERSE_JUMP_TABLE = new int[Byte.MAX_VALUE];
		for (int i = 0; i < CGAT_JUMP_TABLE.length; i++) {
			CGAT_JUMP_TABLE[i] = -1;
			CGAT_REVERSE_JUMP_TABLE[i] = -1;
		}
		CGAT_JUMP_TABLE['C'] = 0;
		CGAT_JUMP_TABLE['G'] = 1;
		CGAT_JUMP_TABLE['A'] = 2;
		CGAT_JUMP_TABLE['T'] = 3;

		CGAT_REVERSE_JUMP_TABLE['C'] = 1;
		CGAT_REVERSE_JUMP_TABLE['G'] = 0;
		CGAT_REVERSE_JUMP_TABLE['A'] = 3;
		CGAT_REVERSE_JUMP_TABLE['T'] = 2;

		DECODE_TABLE = new byte[] { 'C', 'G', 'A', 'T' };
		REVERSE_DECODE_TABLE = new byte[] { 'G', 'C', 'T', 'A' };
	}

	private final int factor;
	private final int len;
	private final int arrSize;
	private transient Object[] root;
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

		this.factor = factor;
		this.len = len;
		this.allowDoubleEntry = allowDoubleEntry;

		int size = 1;
		for (int i = 0; i < factor; i++) {
			size *= 4;
		}
		root = new Object[size];
		arrSize = size;
	}

	private void writeObject(ObjectOutputStream out) throws IOException {
		out.defaultWriteObject();
		writeTree(root, out);
	}

	private void writeTree(Object node, ObjectOutputStream out) throws IOException {
		if (node instanceof Object[]) {
			out.writeByte('[');
			Object[] array = (Object[]) node;
			for (int i = 0; i < arrSize; i++) {
				Object next = array[i];
				if (next != null) {
					out.writeByte('(');
					writeTree(next, out);
				} else {
					out.writeByte('*');
				}
			}
		} else if (node instanceof InternalNullMarker) {
			out.writeByte('x');			
		} else {
			out.writeByte('o');
			out.writeObject(node);
		}
	}

	private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException {
		in.defaultReadObject();
		root = (Object[]) readTree(in);
	}

	private Object readTree(ObjectInputStream in) throws IOException, ClassNotFoundException {
		byte b = in.readByte();
		switch (b) {
		case '[':
			Object[] res = new Object[arrSize];
			for (int i = 0; i < arrSize; i++) {
				res[i] = readTree(in);
			}
			return res;
		case '*':
			return null;
		case 'x':
			return NULL;
		case 'o':
			return in.readObject();
		default:
			throw new IOException("Inconsistent serialization format for kmer trie.");
		}
	}

	public int getLen() {
		return len;
	}

	public long getEntries() {
		return entries;
	}

	public void put(CGATRingBuffer buffer, V value, boolean reverse) {
		if (compressed) {
			throw new IllegalStateException("Cant insert in compressed trie");
		}

		Object[] node = root;

		int mult;
		int pos = 0;
		int j = 0;
		int[] jt = reverse ? CGAT_REVERSE_JUMP_TABLE : CGAT_JUMP_TABLE;

		for (int i = 0; i < len; i += factor) {
			pos = 0;
			mult = 1;
			for (j = 0; j < factor && i + j < len; j++) {
				byte c = buffer.get(reverse ? (len - i - j - 1) : (i + j));
				if (c < 0 || jt[c] == -1) {
					throw new IllegalArgumentException("Not a CGAT sequence");
				}
				pos += jt[c] * mult;
				mult *= 4;
			}
			if (i < len - 1) {
				Object[] next = (Object[]) node[pos];
				if (next == null) {
					next = new Object[root.length];
					node[pos] = next;
				}
				node = next;
			}
		}
		if (!allowDoubleEntry && node[pos] != null && !node[pos].equals(value)) {
			throw new IllegalStateException("double entry " + node[pos] + " " + value);
		}
		if (node[pos] == null) {
			entries++;
		}
		node[pos] = value == null ? NULL : value;
	}

	public void put(byte[] nseq, int start, V value, boolean reverse) {
		if (compressed) {
			throw new IllegalStateException("Cant insert in compressed trie");
		}

		Object[] node = root;

		int mult;
		int pos = 0;
		int j = 0;
		int[] jt = reverse ? CGAT_REVERSE_JUMP_TABLE : CGAT_JUMP_TABLE;

		for (int i = 0; i < len; i += factor) {
			pos = 0;
			mult = 1;
			for (j = 0; j < factor && i + j < len; j++) {
				byte c = CGAT.cgatToUpperCase(nseq[start + (reverse ? (len - i - j - 1) : (i + j))]);
				if (c < 0 || jt[c] == -1) {
					throw new IllegalArgumentException("Not a CGAT sequence");
				}
				pos += jt[c] * mult;
				mult *= 4;
			}
			if (i < len - 1) {
				Object[] next = (Object[]) node[pos];
				if (next == null) {
					next = new Object[root.length];
					node[pos] = next;
				}
				node = next;
			}
		}
		if (!allowDoubleEntry && node[pos] != null && !node[pos].equals(value)) {
			throw new IllegalStateException("double entry " + node[pos] + " " + value);
		}
		if (node[pos] == null) {
			entries++;
		}
		node[pos] = value == null ? NULL : value;
	}

	public void visit(KMerTrieVisitor<V> visitor, boolean reverse) {
		if (compressed) {
			throw new IllegalStateException("Cant collect values on compressed trie (yet)");
		}
		collectValuesHelp(root, 0, new byte[len], visitor, reverse);
	}

	@SuppressWarnings("unchecked")
	private void collectValuesHelp(Object node, int pos, byte[] kmer, KMerTrieVisitor<V> visitor, boolean reverse) {
		if (node == null) {
			return;
		} else if (pos >= len) {
			visitor.nextValue(this, kmer, node instanceof InternalNullMarker ? null : (V) node);
		} else if (node instanceof Object[]) {
			Object[] arr = (Object[]) node;

			for (int i = 0; i < arr.length; i++) {
				int div = 1;
				for (int j = 0; j < factor && pos + j < len; j++) {
					kmer[reverse ? (len - 1 - pos - j)
							: (pos + j)] = (reverse ? REVERSE_DECODE_TABLE : DECODE_TABLE)[(i / div) % 4];
					div *= 4;
				}
				collectValuesHelp(arr[i], pos + factor, kmer, visitor, reverse);
			}
		} else {
			throw new IllegalStateException("no value / no tree node");
		}
	}

	@SuppressWarnings("unchecked")
	public V get(CGATRingBuffer buffer, boolean reverse) {
		Object node = root;

		int mult;
		int pos = 0;
		int j = 0;
		int snPosIndex = 0;
		int[] jt = reverse ? CGAT_REVERSE_JUMP_TABLE : CGAT_JUMP_TABLE;
		int base = (reverse ? len - 1 : 0);

		for (int i = 0; i < len && node != null; i += factor) {
			pos = 0;
			mult = 1;
			for (j = 0; j < factor && i + j < len; j++) {
				byte c = buffer.get(reverse ? (base - i - j) : (i + j));
				if (c < 0 || jt[c] == -1) {
					throw new IllegalArgumentException("Not a CGAT sequence");
				}
				pos += jt[c] * mult;
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
		return node instanceof InternalNullMarker ? null : (V) node;
	}

	@SuppressWarnings("unchecked")
	public V get(byte[] nseq, int start, boolean reverse) {
		Object node = root;

		int mult;
		int pos = 0;
		int j = 0;
		int snPosIndex = 0;
		int[] jt = reverse ? CGAT_REVERSE_JUMP_TABLE : CGAT_JUMP_TABLE;
		int base = start + (reverse ? len - 1 : 0);

		for (int i = 0; i < len && node != null; i += factor) {
			pos = 0;
			mult = 1;
			for (j = 0; j < factor && i + j < len; j++) {
				byte c = CGAT.CGAT_TO_UPPER_CASE[128 + nseq[base + (reverse ? (-i - j) : (i + j))]];
				if (c < 0 || jt[c] == -1) {
					return null; // throw new IllegalArgumentException("Not a CGAT sequence");
				}
				pos += jt[c] * mult;
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
		return node instanceof InternalNullMarker ? null : (V) node;
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
		private transient Object next;

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
			switch (posIndex) {
			case 0:
				return pos0;
			case 1:
				return pos1;
			case 2:
				return pos2;
			case 3:
				return pos3;
			case 4:
				return pos4;
			case 5:
				return pos5;
			case 6:
				return pos6;
			case 7:
				return pos7;
			default:
				throw new IllegalArgumentException("Invalid pos index");
			}
		}

		public int nextPos(int posIndex) {
			int next = (posIndex + 1) % 8;
			return getPosVal(next) == -1 ? 0 : next;
		}
		
		private void writeObject(ObjectOutputStream out) throws IOException {
			out.defaultWriteObject();
			out.writeObject(next);
		}

		private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException {
			in.defaultReadObject();
			next = in.readObject();
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
		ObjectOutputStream oOut = new ObjectOutputStream(StreamProvider.getOutputStreamForFile(trieFile));
		oOut.writeObject(this);
		oOut.close();
	}

	public static KMerTrie<?> load(File filterFile) throws IOException, ClassNotFoundException {
		ObjectInputStream oOut = new ObjectInputStream(StreamProvider.getInputStreamForFile(filterFile));
		KMerTrie<?> res = (KMerTrie<?>) oOut.readObject();
		oOut.close();
		return res;
	}

	public interface KMerTrieVisitor<V extends Serializable> {
		public void nextValue(KMerTrie<V> trie, byte[] kmer, V value);
	}

	private static final class InternalNullMarker implements Serializable {
		private static final long serialVersionUID = 1L;

		private InternalNullMarker() {
		}
	}
}
