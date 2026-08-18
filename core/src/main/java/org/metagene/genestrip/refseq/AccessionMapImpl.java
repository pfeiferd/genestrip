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
 * An {@link AccessionMap} that stores the accession keys in fixed-size slots of a flat arena which,
 * once {@link #optimize()} is called, is sorted by accession key so that lookups use binary search.
 *
 * @deprecated Superseded by {@link AccessionMapTrieImpl}, which holds the same entries in about a
 *             quarter less memory and sorts them about three times faster, by keeping only the part
 *             of a key that its bucket does not already imply. This implementation is kept because
 *             it makes no assumption whatsoever about what an accession looks like: it neither
 *             groups keys by their leading characters nor treats any alphabet as the expected one,
 *             which is what a catalog of a wholly different shape would want.
 */
@Deprecated
public class AccessionMapImpl implements AccessionMap {
	// A key is held in a fixed-size slot of a flat arena rather than in a byte array object of its
	// own: an accession costs the slot alone, where a byte array of ten bytes occupies some thirty
	// once its object header and padding are counted, and filling the map allocates nothing per entry.
	//
	// The slot is what lets the sort work on the keys in place, so no array of positions into the
	// arena is needed at all - which not only saves eight bytes an entry but keeps the sort reading
	// the arena in its own index order, where a permuted array of positions would send every single
	// comparison to a random address.
	//
	// A slot is a length byte followed by the key. A key too long for its slot keeps as much of itself
	// there as fits and puts the whole of itself in the overflow list, whose index occupies the last
	// four bytes of the slot. Two such keys are therefore still ordered by their length and their
	// beginning without consulting the list, which keeps a catalog of consistently long accessions -
	// where every slot would overflow - sorting at nearly the speed of one that fits. Correctness
	// never depends on the slot being wide enough at all. The longest accession in the RefSeq catalog
	// is 20 bytes, so with 23 the overflow list is expected to stay empty.
	//
	// The one length the slot cannot express is one that does not fit its length byte, so a key of
	// more than 255 bytes is rejected outright rather than stored in a way the order could not
	// reproduce. No accession comes anywhere near that.
	private static final int MAX_INLINE_KEY_LENGTH = 23;
	private static final int STRIDE = 1 + MAX_INLINE_KEY_LENGTH;
	private static final int OVERFLOW_INDEX_BYTES = Integer.BYTES;
	private static final int OVERFLOW_INDEX_OFFSET = STRIDE - OVERFLOW_INDEX_BYTES;
	private static final int OVERFLOW_PREFIX_LENGTH = OVERFLOW_INDEX_OFFSET - 1;
	private static final int MAX_KEY_LENGTH = 255;
	// A slot is read, written and compared eight bytes at a time: the stride is a multiple of eight,
	// so a slot is exactly three long words. Reading them big-endian makes the unsigned order of the
	// words the order of the bytes, and since a slot begins with the length byte and every byte after
	// the key is zero - a slot is written once and the sort only ever moves whole slots - comparing
	// the words of two slots yields exactly the intended order: shorter keys first, then by content.
	private static final VarHandle LONG_VIEW =
			MethodHandles.byteArrayViewVarHandle(long[].class, ByteOrder.BIG_ENDIAN);

	// The arena and the values live in chunks that are added as they fill up, rather than in one array
	// that has to be reallocated and copied. Filling therefore needs no separate counting pass over
	// the (multi-gigabyte) catalog to size the map, without a growing copy ever holding the old and
	// the new storage at once. Chunk sizes start at about 16 MB and double up to about 256 MB, so a
	// map with few entries stays small while a large one still ends up in a handful of chunks - which
	// keeps both the unused tail of the last chunk and the per-chunk overhead of the garbage collector
	// negligible. The slot count of a chunk is a power of two, so an entry index splits into chunk and
	// slot by a shift and a mask beyond the doubling ramp, and by a base-two logarithm inside it.
	private static final int BASE_BITS = 19;
	private static final int MAX_BITS = 23;
	private static final int STEPS = MAX_BITS - BASE_BITS;
	private static final int RAMP_END = (1 << BASE_BITS) * ((1 << STEPS) - 1);
	private static final int MAX_SIZE = 1 << MAX_BITS;
	private static final int MAX_MASK = MAX_SIZE - 1;
	// The chunk index itself is a plain array rather than a list: it holds one reference per chunk, so
	// even a few dozen slots cover far more entries than an int index can address. It is doubled if
	// that should ever not be enough.
	private static final int INITIAL_CHUNKS = 100;

	private byte[][] arena;
	private TaxIdNode[][] valueChunks;
	private int chunks;
	private int chunkStart;
	private int capacity;

	private List<byte[]> overflow;

	private int entries;
	private boolean sorted;

	/**
	 * Creates a map that takes as many entries as are put into it.
	 */
	public AccessionMapImpl() {
		entries = 0;
		arena = new byte[INITIAL_CHUNKS][];
		valueChunks = new TaxIdNode[INITIAL_CHUNKS][];
		chunks = 0;
		chunkStart = 0;
		capacity = 0;
		overflow = null;
		sorted = false;
	}

	/**
	 * Creates a map, reserving room for roughly the given number of entries. The size is a hint only -
	 * the map takes as many entries as are put into it either way.
	 *
	 * @param size the number of entries to reserve room for
	 */
	public AccessionMapImpl(int size) {
		this();
		while (capacity < size) {
			addChunk();
		}
	}

	/**
	 * Returns the number of slots the chunk with the given index holds.
	 */
	private static int chunkSize(int chunk) {
		return chunk < STEPS ? 1 << (BASE_BITS + chunk) : MAX_SIZE;
	}

	/**
	 * Splits an entry index into the chunk holding it and the slot within that chunk, returned as
	 * {@code (chunk << 32) | slot}.
	 */
	private static long locate(int index) {
		if (index >= RAMP_END) {
			final int beyond = index - RAMP_END;
			return (((long) STEPS + (beyond >>> MAX_BITS)) << 32) | (beyond & MAX_MASK);
		}
		final int chunk = 31 - Integer.numberOfLeadingZeros((index >>> BASE_BITS) + 1);
		final int start = (1 << BASE_BITS) * ((1 << chunk) - 1);
		return ((long) chunk << 32) | (index - start);
	}

	private void addChunk() {
		if (chunks == arena.length) {
			arena = Arrays.copyOf(arena, arena.length * 2);
			valueChunks = Arrays.copyOf(valueChunks, valueChunks.length * 2);
		}
		final int size = chunkSize(chunks);
		arena[chunks] = new byte[size * STRIDE];
		valueChunks[chunks] = new TaxIdNode[size];
		chunkStart = capacity;
		capacity += size;
		chunks++;
	}

	private TaxIdNode valueAt(int index) {
		final long location = locate(index);
		return valueChunks[(int) (location >>> 32)][(int) location];
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
		int pos = binaryKeySearch(0, entries, array, start, end);
		if (pos < 0) {
			return null;
		}
		return valueAt(pos);
	}

	public int getEntriesForNode(TaxIdNode node) {
		int counter = 0;
		for (int i = 0; i < entries; i++) {
			if (valueAt(i) == node) {
				counter++;
			}
		}
		return counter;
	}

	/**
	 * Binary-searches the sorted keys in {@code [from, to)} for the key {@code array[start, end)},
	 * returning its index, or {@code -(insertionPoint + 1)} if absent.
	 *
	 * @param from  the index of the first key to search (inclusive)
	 * @param to    the index one past the last key to search (exclusive)
	 * @param array the byte array holding the search key
	 * @param start the start index of the search key in {@code array} (inclusive)
	 * @param end   the end index of the search key in {@code array} (exclusive)
	 * @return the index of the matching key, or {@code -(insertionPoint + 1)} if absent
	 */
	protected int binaryKeySearch(int from, int to, byte[] array, int start, int end) {
		to--;
		while (from <= to) {
			final int mid = (from + to) >>> 1;
			final long location = locate(mid);
			final byte[] chunk = arena[(int) (location >>> 32)];
			final int offset = (int) location * STRIDE;
			final int res = compareKey(chunk, offset, array, start, end - start);
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
		if (entries == capacity) {
			addChunk();
		}
		final int slot = entries - chunkStart;
		final byte[] chunk = arena[chunks - 1];
		final int offset = slot * STRIDE;
		final int len = end - start;
		if (len <= MAX_INLINE_KEY_LENGTH) {
			chunk[offset] = (byte) len;
			System.arraycopy(array, start, chunk, offset + 1, len);
		} else if (len > MAX_KEY_LENGTH) {
			throw new IllegalArgumentException("Accession key of length " + len
					+ " exceeds the maximum " + MAX_KEY_LENGTH + " this map can store.");
		} else {
			// Too long for a slot: the whole key goes to the overflow list, while the slot keeps its
			// length, as much of its beginning as fits and the index of the key in that list.
			if (overflow == null) {
				overflow = new ArrayList<>();
			}
			final int index = overflow.size();
			overflow.add(Arrays.copyOfRange(array, start, end));
			chunk[offset] = (byte) len;
			System.arraycopy(array, start, chunk, offset + 1, OVERFLOW_PREFIX_LENGTH);
			chunk[offset + OVERFLOW_INDEX_OFFSET] = (byte) (index >>> 24);
			chunk[offset + OVERFLOW_INDEX_OFFSET + 1] = (byte) (index >>> 16);
			chunk[offset + OVERFLOW_INDEX_OFFSET + 2] = (byte) (index >>> 8);
			chunk[offset + OVERFLOW_INDEX_OFFSET + 3] = (byte) index;
		}
		valueChunks[chunks - 1][slot] = node;
		entries++;
	}

	public void optimize() {
		// The keys are sorted where they lie, with the values permuted in lock-step, so that a lookup
		// can binary-search them afterwards.
		it.unimi.dsi.fastutil.Arrays.quickSort(0, entries, new IntComparator() {
			@Override
			public int compare(int k1, int k2) {
				final long location1 = locate(k1);
				final long location2 = locate(k2);
				return compareKeys(arena[(int) (location1 >>> 32)], (int) location1 * STRIDE,
						arena[(int) (location2 >>> 32)], (int) location2 * STRIDE);
			}
		}, new Swapper() {
			@Override
			public void swap(int a, int b) {
				final long locationA = locate(a);
				final long locationB = locate(b);
				final int chunkA = (int) (locationA >>> 32);
				final int chunkB = (int) (locationB >>> 32);
				final int slotA = (int) locationA;
				final int slotB = (int) locationB;
				final byte[] keysA = arena[chunkA];
				final byte[] keysB = arena[chunkB];
				final int offsetA = slotA * STRIDE;
				final int offsetB = slotB * STRIDE;
				for (int i = 0; i < STRIDE; i += Long.BYTES) {
					final long word = (long) LONG_VIEW.get(keysA, offsetA + i);
					LONG_VIEW.set(keysA, offsetA + i, (long) LONG_VIEW.get(keysB, offsetB + i));
					LONG_VIEW.set(keysB, offsetB + i, word);
				}
				final TaxIdNode value = valueChunks[chunkA][slotA];
				valueChunks[chunkA][slotA] = valueChunks[chunkB][slotB];
				valueChunks[chunkB][slotB] = value;
			}
		});
		sorted = true;
	}

	/**
	 * Returns the key held in the overflow list for the slot at the given offset.
	 */
	private byte[] overflowKey(byte[] chunk, int offset) {
		final int index = ((chunk[offset + OVERFLOW_INDEX_OFFSET] & 0xFF) << 24)
				| ((chunk[offset + OVERFLOW_INDEX_OFFSET + 1] & 0xFF) << 16)
				| ((chunk[offset + OVERFLOW_INDEX_OFFSET + 2] & 0xFF) << 8)
				| (chunk[offset + OVERFLOW_INDEX_OFFSET + 3] & 0xFF);
		return overflow.get(index);
	}

	/**
	 * Lexicographically compares the keys in the two given slots, ordering shorter keys before longer
	 * ones.
	 *
	 * @param chunk1  the arena chunk holding the first slot
	 * @param offset1 the offset of the first slot in its chunk
	 * @param chunk2  the arena chunk holding the second slot
	 * @param offset2 the offset of the second slot in its chunk
	 * @return a negative, zero or positive value if the first key sorts before, equal to or after the second
	 */
	protected int compareKeys(byte[] chunk1, int offset1, byte[] chunk2, int offset2) {
		final int len1 = chunk1[offset1] & 0xFF;
		final int len2 = chunk2[offset2] & 0xFF;
		if (len1 <= MAX_INLINE_KEY_LENGTH && len2 <= MAX_INLINE_KEY_LENGTH) {
			for (int i = 0; i < STRIDE; i += Long.BYTES) {
				final long word1 = (long) LONG_VIEW.get(chunk1, offset1 + i);
				final long word2 = (long) LONG_VIEW.get(chunk2, offset2 + i);
				if (word1 != word2) {
					return Long.compareUnsigned(word1, word2);
				}
			}
			return 0;
		}
		// At least one key is too long for its slot. The stored lengths order the two, and where those
		// agree the stored beginnings do, so the overflow list is reached only when both agree.
		if (len1 != len2) {
			return len1 - len2;
		}
		final int prefix = compareBytes(chunk1, offset1 + 1, OVERFLOW_PREFIX_LENGTH,
				chunk2, offset2 + 1, OVERFLOW_PREFIX_LENGTH);
		if (prefix != 0) {
			return prefix;
		}
		final byte[] key1 = overflowKey(chunk1, offset1);
		final byte[] key2 = overflowKey(chunk2, offset2);
		return compareBytes(key1, 0, key1.length, key2, 0, key2.length);
	}

	/**
	 * Lexicographically compares the key in the given slot against {@code array[start, start + len)},
	 * ordering the shorter of the two before the longer.
	 *
	 * @param chunk  the arena chunk holding the slot
	 * @param offset the offset of the slot in its chunk
	 * @param array  the byte array holding the other key
	 * @param start  the start index of the other key in {@code array} (inclusive)
	 * @param len    the length of the other key
	 * @return a negative, zero or positive value if the key sorts before, equal to or after the other key
	 */
	protected int compareKey(byte[] chunk, int offset, byte[] array, int start, int len) {
		final int keyLen = chunk[offset] & 0xFF;
		if (keyLen <= MAX_INLINE_KEY_LENGTH) {
			return compareBytes(chunk, offset + 1, keyLen, array, start, len);
		}
		final byte[] key = overflowKey(chunk, offset);
		return compareBytes(key, 0, key.length, array, start, len);
	}

	/**
	 * Lexicographically compares two byte ranges, ordering the shorter of the two before the longer.
	 *
	 * @param array1 the byte array holding the first key
	 * @param start1 the start index of the first key (inclusive)
	 * @param len1   the length of the first key
	 * @param array2 the byte array holding the second key
	 * @param start2 the start index of the second key (inclusive)
	 * @param len2   the length of the second key
	 * @return a negative, zero or positive value if the first key sorts before, equal to or after the second
	 */
	protected static int compareBytes(byte[] array1, int start1, int len1, byte[] array2, int start2, int len2) {
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
}
