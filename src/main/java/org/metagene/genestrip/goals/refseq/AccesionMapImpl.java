package org.metagene.genestrip.goals.refseq;

import java.util.Arrays;

import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

import it.unimi.dsi.fastutil.BigArrays;
import it.unimi.dsi.fastutil.BigSwapper;
import it.unimi.dsi.fastutil.longs.LongComparator;

public class AccesionMapImpl implements AccessionMap {
	private final byte[][] keys;
	private final TaxIdNode[] values;

	private int entries;
	private boolean sorted;

	public AccesionMapImpl(int size) {
		entries = 0;
		keys = new byte[size][];
		values = new TaxIdNode[size];
		sorted = false;
	}

	@Override
	public TaxIdNode get(byte[] array, int start, int end) {
		if (!sorted) {
			throw new IllegalStateException("Map must be optimized before get.");
		}
		int pos = binaryKeySearch(0, entries, array, start, end);
		if (pos < 0) {
			return null;
		}
		return values[pos];
	}

	protected int binaryKeySearch(int from, int to, byte[] array, int start, int end) {
		byte[] midVal;
		to--;
		while (from <= to) {
			final int mid = (from + to) >>> 1;
			midVal = keys[mid];
			final int res = compareBytes(midVal, array, start, end);
			if (res < 0)
				from = mid + 1;
			else if (res > 0)
				to = mid - 1;
			else
				return mid;
		}
		return -(from + 1);

	}

	@Override
	public void put(byte[] array, int start, int end, TaxIdNode node) {
		sorted = false;
		byte[] key = Arrays.copyOfRange(array, start, end);
		keys[entries] = key;
		values[entries] = node;
		entries++;
	}

	public void optimize() {
		BigArrays.quickSort(0, entries, new LongComparator() {
			@Override
			public int compare(long k1, long k2) {
				byte[] key1 = keys[(int) k1];
				byte[] key2 = keys[(int) k2];
				return compareBytes(key1, key2);
			}
		}, new BigSwapper() {
			@Override
			public void swap(long a, long b) {
				byte[] key1 = keys[(int) a];
				byte[] key2 = keys[(int) b];
				keys[(int) b] = key1;
				keys[(int) a] = key2;
				TaxIdNode n1 = values[(int) a];
				TaxIdNode n2 = values[(int) b];
				values[(int) b] = n1;
				values[(int) a] = n2;
			}
		});
		sorted = true;
	}

	protected int compareBytes(byte[] key1, byte[] key2) {
		if (key1.length == key2.length) {
			for (int i = 0; i < key1.length; i++) {
				if (key1[i] != key2[i]) {
					return key1[i] - key2[i];
				}
			}
			return 0;
		} else {
			return key1.length - key2.length;
		}
	}

	protected int compareBytes(byte[] key1, byte[] array, int start, int end) {
		int len = end - start;
		if (key1.length == len) {
			for (int i = 0; i < len; i++) {
				if (key1[i] != array[start + i]) {
					return key1[i] - array[start + i];
				}
			}
			return 0;
		} else {
			return key1.length - len;
		}
	}
}
