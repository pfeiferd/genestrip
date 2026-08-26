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
package org.metagene.genestrip.refseq;

import java.lang.invoke.MethodHandles;
import java.lang.invoke.VarHandle;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

import it.unimi.dsi.fastutil.Swapper;
import it.unimi.dsi.fastutil.ints.IntComparator;

/**
 * An {@link AccessionMap} that dispatches on the first characters of an accession through a trie of
 * jump tables and keeps the remainder of the key in a fixed-size slot of a flat arena per bucket,
 * which {@link #optimize()} sorts so that lookups binary-search a single bucket.
 * <p>
 * Compared to {@link AccessionMapImpl} this trades a handful of table lookups for a shorter key and
 * a much shorter binary search: the leading characters are implied by the bucket a key sits in and
 * are not stored, and a bucket holds a small enough share of the entries that its arena is
 * comparable to the size of a processor cache.
 */
public class AccessionMapTrieImpl implements AccessionMap {
	// The first characters of a key select a bucket by walking one jump table per character. A table
	// per character is what keeps the trie small: a single flat table over the same characters would
	// need one slot per possible combination, of which all but a few thousand stay empty.
	private static final int PREFIX_LENGTH = 6;
	// The characters an accession is made of - digits, letters and the two punctuation marks - which
	// is what the jump tables are sized for. Anything else is not wrong, only slower: such a key goes
	// to the fallback bucket, which stores keys whole.
	private static final String ALPHABET = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
			+ "abcdefghijklmnopqrstuvwxyz._";
	private static final int ALPHABET_SIZE = ALPHABET.length();
	private static final byte[] CHARACTER_INDEX = new byte[256];

	// A slot holds the length of the stored key remainder followed by the remainder itself. As in
	// AccessionMapImpl, a remainder too long for its slot keeps as much of itself there as fits and
	// the index of the whole of it in the bucket's overflow list in the last four bytes, and the
	// length byte cannot express more than 255.
	private static final int STRIDE = 16;
	private static final int MAX_INLINE_TAIL_LENGTH = STRIDE - 1;
	private static final int MAX_TAIL_LENGTH = 255;
	private static final int OVERFLOW_INDEX_OFFSET = STRIDE - Integer.BYTES;
	private static final int OVERFLOW_PREFIX_LENGTH = OVERFLOW_INDEX_OFFSET - 1;
	// A slot is exactly two long words, read big-endian so that the unsigned order of the words is
	// the order of the bytes - the length byte leading, so shorter remainders sort first.
	private static final VarHandle LONG_VIEW =
			MethodHandles.byteArrayViewVarHandle(long[].class, ByteOrder.BIG_ENDIAN);

	private static final int INITIAL_BUCKET_SLOTS = 16;

	static {
		Arrays.fill(CHARACTER_INDEX, (byte) -1);
		for (int i = 0; i < ALPHABET_SIZE; i++) {
			CHARACTER_INDEX[ALPHABET.charAt(i)] = (byte) i;
		}
	}

	private final Object[] root;
	private final List<Bucket> buckets;
	private Bucket fallback;

	private int entries;
	private boolean sorted;

	/**
	 * Creates a map that takes as many entries as are put into it.
	 */
	/** The genomes a database built from this map is to hold, or null if no genome limit is in force. */
	private final GenomeKeyTrie admittedGenomes;

	/**
	 * Creates an empty map holding no genome selection, for a database that is to hold every genome.
	 */
	public AccessionMapTrieImpl() {
		this(null);
	}

	/**
	 * Creates an empty map carrying the given genome selection.
	 *
	 * @param admittedGenomes the genomes a database built from this map is to hold, or {@code null}
	 *                        where {@code maxGenomesPerTaxid} is not set
	 */
	public AccessionMapTrieImpl(GenomeKeyTrie admittedGenomes) {
		this.admittedGenomes = admittedGenomes;
		root = new Object[ALPHABET_SIZE];
		buckets = new ArrayList<>();
		fallback = null;
		entries = 0;
		sorted = false;
	}

	/**
	 * The entries of one prefix, holding the remainders of their keys and their values.
	 */
	private static class Bucket {
		private final int stripped;
		private byte[] arena;
		private TaxIdNode[] values;
		private int size;
		private List<byte[]> overflow;

		Bucket(int stripped) {
			this.stripped = stripped;
			arena = new byte[INITIAL_BUCKET_SLOTS * STRIDE];
			values = new TaxIdNode[INITIAL_BUCKET_SLOTS];
			size = 0;
		}

		/**
		 * Makes room for one more entry, growing by half so that the copying stays proportional to
		 * the bucket rather than to the map, and never holds two copies of more than one bucket.
		 */
		private void ensureCapacity() {
			if (size == values.length) {
				final int capacity = values.length + (values.length >> 1) + 1;
				arena = Arrays.copyOf(arena, capacity * STRIDE);
				values = Arrays.copyOf(values, capacity);
			}
		}

		/**
		 * Releases the room the bucket reserved but did not use.
		 */
		private void trim() {
			if (size != values.length) {
				arena = Arrays.copyOf(arena, size * STRIDE);
				values = Arrays.copyOf(values, size);
			}
		}
	}

	/**
	 * Returns the bucket the given key belongs to, creating it if asked, or {@code null} if it does
	 * not exist. Keys shorter than the prefix, or holding a character the jump tables do not cover,
	 * belong to the fallback bucket, which strips nothing.
	 */
	private Bucket bucketFor(byte[] array, int start, int end, boolean create) {
		if (end - start < PREFIX_LENGTH) {
			return fallbackBucket(create);
		}
		Object[] node = root;
		for (int i = 0; i < PREFIX_LENGTH; i++) {
			final int index = CHARACTER_INDEX[array[start + i] & 0xFF];
			if (index < 0) {
				return fallbackBucket(create);
			}
			if (i == PREFIX_LENGTH - 1) {
				Bucket bucket = (Bucket) node[index];
				if (bucket == null && create) {
					bucket = new Bucket(PREFIX_LENGTH);
					node[index] = bucket;
					buckets.add(bucket);
				}
				return bucket;
			}
			Object[] child = (Object[]) node[index];
			if (child == null) {
				if (!create) {
					return null;
				}
				child = new Object[ALPHABET_SIZE];
				node[index] = child;
			}
			node = child;
		}
		throw new IllegalStateException();
	}

	private Bucket fallbackBucket(boolean create) {
		if (fallback == null && create) {
			fallback = new Bucket(0);
			buckets.add(fallback);
		}
		return fallback;
	}

	@Override
	public void put(byte[] array, int start, int end, TaxIdNode node) {
		sorted = false;
		final Bucket bucket = bucketFor(array, start, end, true);
		final int from = start + bucket.stripped;
		final int len = end - from;
		if (len > MAX_TAIL_LENGTH) {
			throw new IllegalArgumentException("Accession key of length " + (end - start)
					+ " exceeds the maximum " + (MAX_TAIL_LENGTH + bucket.stripped)
					+ " this map can store.");
		}
		bucket.ensureCapacity();
		final int offset = bucket.size * STRIDE;
		final byte[] arena = bucket.arena;
		if (len <= MAX_INLINE_TAIL_LENGTH) {
			arena[offset] = (byte) len;
			System.arraycopy(array, from, arena, offset + 1, len);
		} else {
			if (bucket.overflow == null) {
				bucket.overflow = new ArrayList<>();
			}
			final int index = bucket.overflow.size();
			bucket.overflow.add(Arrays.copyOfRange(array, from, end));
			arena[offset] = (byte) len;
			System.arraycopy(array, from, arena, offset + 1, OVERFLOW_PREFIX_LENGTH);
			arena[offset + OVERFLOW_INDEX_OFFSET] = (byte) (index >>> 24);
			arena[offset + OVERFLOW_INDEX_OFFSET + 1] = (byte) (index >>> 16);
			arena[offset + OVERFLOW_INDEX_OFFSET + 2] = (byte) (index >>> 8);
			arena[offset + OVERFLOW_INDEX_OFFSET + 3] = (byte) index;
		}
		bucket.values[bucket.size] = node;
		bucket.size++;
		entries++;
	}

	public void optimize() {
		// Every bucket is sorted on its own: keys in one bucket share their leading characters, so
		// ordering their remainders orders the keys themselves, and no order holds between buckets
		// because a lookup never sees more than one of them.
		for (Bucket bucket : buckets) {
			bucket.trim();
			sort(bucket);
		}
		sorted = true;
	}

	private void sort(final Bucket bucket) {
		final byte[] arena = bucket.arena;
		final TaxIdNode[] values = bucket.values;
		it.unimi.dsi.fastutil.Arrays.quickSort(0, bucket.size, new IntComparator() {
			@Override
			public int compare(int k1, int k2) {
				return compareSlots(bucket, arena, k1 * STRIDE, arena, k2 * STRIDE);
			}
		}, new Swapper() {
			@Override
			public void swap(int a, int b) {
				final int offsetA = a * STRIDE;
				final int offsetB = b * STRIDE;
				for (int i = 0; i < STRIDE; i += Long.BYTES) {
					final long word = (long) LONG_VIEW.get(arena, offsetA + i);
					LONG_VIEW.set(arena, offsetA + i, (long) LONG_VIEW.get(arena, offsetB + i));
					LONG_VIEW.set(arena, offsetB + i, word);
				}
				final TaxIdNode value = values[a];
				values[a] = values[b];
				values[b] = value;
			}
		});
	}

	@Override
	public TaxIdNode get(byte[] array, int start, int end, boolean assemblyAccessionsOnly) {
		if (!sorted) {
			throw new IllegalStateException("Map must be optimized before get.");
		}
		if (assemblyAccessionsOnly &&
				!AccessionFileProcessor.isAssemblyAccession(array, start) &&
				!AccessionFileProcessor.isMRNAAccession(array, start) &&
				!AccessionFileProcessor.isRNAAccession(array, start)) {
			return null;
		}
		final Bucket bucket = bucketFor(array, start, end, false);
		if (bucket == null) {
			return null;
		}
		final int from = start + bucket.stripped;
		int low = 0;
		int high = bucket.size - 1;
		while (low <= high) {
			final int mid = (low + high) >>> 1;
			final int res = compareSlotToKey(bucket, bucket.arena, mid * STRIDE, array, from, end - from);
			if (res < 0)
				low = mid + 1;
			else if (res > 0)
				high = mid - 1;
			else
				return bucket.values[mid];
		}
		return null;
	}

	public int getEntriesForNode(TaxIdNode node) {
		int counter = 0;
		for (Bucket bucket : buckets) {
			for (int i = 0; i < bucket.size; i++) {
				if (bucket.values[i] == node) {
					counter++;
				}
			}
		}
		return counter;
	}

	/**
	 * Returns the key remainder held in the bucket's overflow list for the slot at the given offset.
	 */
	private static byte[] overflowTail(Bucket bucket, byte[] arena, int offset) {
		final int index = ((arena[offset + OVERFLOW_INDEX_OFFSET] & 0xFF) << 24)
				| ((arena[offset + OVERFLOW_INDEX_OFFSET + 1] & 0xFF) << 16)
				| ((arena[offset + OVERFLOW_INDEX_OFFSET + 2] & 0xFF) << 8)
				| (arena[offset + OVERFLOW_INDEX_OFFSET + 3] & 0xFF);
		return bucket.overflow.get(index);
	}

	/**
	 * Lexicographically compares the key remainders in two slots of the same bucket, ordering
	 * shorter ones before longer ones.
	 */
	private static int compareSlots(Bucket bucket, byte[] arena1, int offset1, byte[] arena2, int offset2) {
		final int len1 = arena1[offset1] & 0xFF;
		final int len2 = arena2[offset2] & 0xFF;
		if (len1 <= MAX_INLINE_TAIL_LENGTH && len2 <= MAX_INLINE_TAIL_LENGTH) {
			for (int i = 0; i < STRIDE; i += Long.BYTES) {
				final long word1 = (long) LONG_VIEW.get(arena1, offset1 + i);
				final long word2 = (long) LONG_VIEW.get(arena2, offset2 + i);
				if (word1 != word2) {
					return Long.compareUnsigned(word1, word2);
				}
			}
			return 0;
		}
		if (len1 != len2) {
			return len1 - len2;
		}
		final int prefix = compareBytes(arena1, offset1 + 1, OVERFLOW_PREFIX_LENGTH,
				arena2, offset2 + 1, OVERFLOW_PREFIX_LENGTH);
		if (prefix != 0) {
			return prefix;
		}
		final byte[] tail1 = overflowTail(bucket, arena1, offset1);
		final byte[] tail2 = overflowTail(bucket, arena2, offset2);
		return compareBytes(tail1, 0, tail1.length, tail2, 0, tail2.length);
	}

	/**
	 * Lexicographically compares the key remainder in a slot against {@code array[start, start + len)},
	 * ordering the shorter of the two before the longer.
	 */
	private static int compareSlotToKey(Bucket bucket, byte[] arena, int offset, byte[] array, int start, int len) {
		final int tailLen = arena[offset] & 0xFF;
		if (tailLen <= MAX_INLINE_TAIL_LENGTH) {
			return compareBytes(arena, offset + 1, tailLen, array, start, len);
		}
		final byte[] tail = overflowTail(bucket, arena, offset);
		return compareBytes(tail, 0, tail.length, array, start, len);
	}

	/**
	 * Lexicographically compares two byte ranges, ordering the shorter of the two before the longer.
	 */
	private static int compareBytes(byte[] array1, int start1, int len1, byte[] array2, int start2, int len2) {
		if (len1 != len2) {
			return len1 - len2;
		}
		for (int i = 0; i < len1; i++) {
			final int b1 = array1[start1 + i] & 0xFF;
			final int b2 = array2[start2 + i] & 0xFF;
			if (b1 != b2) {
				return b1 - b2;
			}
		}
		return 0;
	}

	@Override
	public GenomeKeyTrie getAdmittedGenomes() {
		return admittedGenomes;
	}
}
