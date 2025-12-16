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
package org.metagene.genestrip.bloom;

import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.LargeBitVector;

import java.io.*;
import java.util.Random;

public abstract class AbstractKMerBloomFilter implements Serializable {
	public static long MAX_SMALL_CAPACITY = Integer.MAX_VALUE - 8;

	private static final long serialVersionUID = 1L;

	protected final int k;
	protected final double fpp;
	protected final Random random;

	protected long expectedInsertions;
	protected LargeBitVector bitVector;
	protected int hashes;
	protected long[] hashFactors;
	protected long entries;

	public AbstractKMerBloomFilter(int k, double fpp) {
		if (k <= 0) {
			throw new IllegalArgumentException("k-mer length k must be > 0");
		}
		if (fpp <= 0 || fpp >= 1) {
			throw new IllegalArgumentException("fpp must be a probability");
		}
		this.fpp = fpp;
		this.k = k;

		random = new Random(42);

		bitVector = new LargeBitVector(0);
		entries = 0;
	}

	public void clear() {
		bitVector.clear();
		entries = 0;
	}

	public long ensureExpectedSize(long expectedInsertions, boolean enforceLarge) {
		if (expectedInsertions < 0) {
			throw new IllegalArgumentException("expected insertions must be > 0");
		}
		this.expectedInsertions = expectedInsertions;

		long bits = optimalNumOfBits(expectedInsertions, fpp);
		if (bitVector.ensureCapacity(bits, enforceLarge)) {
			hashes = optimalNumOfHashFunctions(expectedInsertions, bits);
			hashFactors = new long[hashes];
			for (int i = 0; i < hashFactors.length; i++) {
				hashFactors[i] = random.nextLong();
			}
		}
		return bits;
	}

	public long getExpectedInsertions() {
		return expectedInsertions;
	}

	public int getK() {
		return k;
	}

	public int getHashes() {
		return hashes;
	}

	public double getFpp() {
		return fpp;
	}

	public long getBitSize() {
		return bitVector.getBitSize();
	}

	public boolean isLarge() {
		return bitVector.isLarge();
	}

	public long getEntries() {
		return entries;
	}

	protected int optimalNumOfHashFunctions(long n, long m) {
		return Math.max(1, (int) Math.round(((double) m) / n * Math.log(2)));
	}

	protected long optimalNumOfBits(long n, double p) {
		return (long) (-n * Math.log(p) / (Math.log(2) * Math.log(2)));
	}

	/**
	 * Competing calls to putLong are multi-threading enabled.
	 * @param data
	 */
	public void putLong(long data) {
		putViaHash(data);
	}

	protected void putViaHash(long data) {
		synchronized (this) {
			entries++;
		}
		for (int i = 0; i < hashes; i++) {
			bitVector.set(hash(data, i));
		}
	}

	public final boolean contains(byte[] seq, int start, int[] badPos) {
		long data = CGAT.kMerToLong(seq, start, k, badPos);
		if (data == -1 && badPos != null && badPos[0] == -1) {
			return false;
		}
		return containsViaHash(data);
	}

	protected final boolean containsViaHash(final long data) {
		for (int i = 0; i < hashes; i++) {
			if (!bitVector.get(hash(data, i))) {
				return false;
			}
		}
		return true;
	}

	public boolean containsLong(final long data) {
		return containsViaHash(data);
	}

	protected abstract long hash(final long data, final int i);

	public void save(File filterFile) throws IOException {
		try (ObjectOutputStream oOut = new ObjectOutputStream(StreamProvider.getOutputStreamForFile(filterFile))) {
			oOut.writeObject(this);			
		}
	}

	public static AbstractKMerBloomFilter load(InputStream is) throws IOException, ClassNotFoundException {
		try (ObjectInputStream oOut = new ObjectInputStream(is)) {
			return (AbstractKMerBloomFilter) oOut.readObject();
		}
	}

	public static AbstractKMerBloomFilter load(File filterFile) throws IOException, ClassNotFoundException {
		try (InputStream is = StreamProvider.getInputStreamForFile(filterFile)) {
			return load(is);
		}
	}
}