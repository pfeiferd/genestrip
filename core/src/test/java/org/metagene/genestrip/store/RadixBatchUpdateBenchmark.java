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
package org.metagene.genestrip.store;

import java.util.Random;

import org.metagene.genestrip.store.KMerStore.UpdateValueProvider;

/**
 * Stand-alone micro-benchmark and correctness check for {@link RadixKMerStore#updateBatch}.
 *
 * It builds a radix store large enough to exceed the CPU cache, then compares the throughput of the
 * per-k-mer {@link RadixKMerStore#update} against the batched {@link RadixKMerStore#updateBatch}
 * (which overlaps independent cache-missing loads across a batch). It also asserts that the batched
 * path yields exactly the same store contents as the sequential path.
 *
 * Run (needs no network) with e.g.:
 *   java -Xmx12g -cp core/target/classes:core/target/test-classes:base/target/classes:<deps> \
 *       org.metagene.genestrip.store.RadixBatchUpdateBenchmark 20000000 20000000 5 1024
 * Args (all optional): entries, queries, reps, batchSize.
 */
public class RadixBatchUpdateBenchmark {
	private static final int K = 31;
	private static final int RADIX_BITS = 17;
	private static final int VALUE_POOL = 1000;
	private static final long KMER_MASK = (1L << (2 * K)) - 1;

	public static void main(String[] args) {
		int entries = args.length > 0 ? Integer.parseInt(args[0]) : 20_000_000;
		int queries = args.length > 1 ? Integer.parseInt(args[1]) : 20_000_000;
		int reps = args.length > 2 ? Integer.parseInt(args[2]) : 5;
		int batchSize = args.length > 3 ? Integer.parseInt(args[3]) : 1024;

		System.out.printf("entries=%d queries=%d reps=%d batchSize=%d%n", entries, queries, reps, batchSize);

		correctnessCheck();

		long[] kmers = distinctKmers(entries, 42);
		RadixKMerStore<String> store = buildStore(kmers);
		System.out.printf("store built: entries=%d nValues=%d%n", store.getEntries(), store.getNValues());

		// 90% hits (recurring stored k-mers, as in a RefSeq re-read), 10% misses.
		long[] queryKmers = buildQueries(kmers, queries, 111);

		UpdateValueProvider<String> noMove = old -> old;

		// Warm up the JIT on both paths.
		runSequential(store, queryKmers, noMove);
		runBatched(store, queryKmers, noMove, batchSize);

		long seqNs = 0, batNs = 0;
		for (int r = 0; r < reps; r++) {
			seqNs += runSequential(store, queryKmers, noMove);
			batNs += runBatched(store, queryKmers, noMove, batchSize);
		}
		double seqPer = (double) seqNs / ((long) reps * queries);
		double batPer = (double) batNs / ((long) reps * queries);
		System.out.printf("%nsequential: %.1f ns/update  (%.1f M updates/s)%n", seqPer, 1000.0 / seqPer);
		System.out.printf("batched   : %.1f ns/update  (%.1f M updates/s)%n", batPer, 1000.0 / batPer);
		System.out.printf("speedup   : %.2fx%n", seqPer / batPer);
	}

	private static long runSequential(RadixKMerStore<String> store, long[] q, UpdateValueProvider<String> p) {
		long t0 = System.nanoTime();
		for (int i = 0; i < q.length; i++) {
			store.update(q[i], p);
		}
		return System.nanoTime() - t0;
	}

	private static long runBatched(RadixKMerStore<String> store, long[] q, UpdateValueProvider<String> p, int batchSize) {
		RadixKMerStore.BatchBuffers buf = new RadixKMerStore.BatchBuffers(batchSize);
		long t0 = System.nanoTime();
		for (int i = 0; i < q.length; i++) {
			if (buf.add(q[i])) {
				store.updateBatch(buf, p);
			}
		}
		if (!buf.isEmpty()) {
			store.updateBatch(buf, p);
		}
		return System.nanoTime() - t0;
	}

	private static void correctnessCheck() {
		int n = 1_000_000;
		long[] kmers = distinctKmers(n, 7);
		RadixKMerStore<String> a = buildStore(kmers);
		RadixKMerStore<String> b = buildStore(kmers);
		long[] q = buildQueries(kmers, 2_000_000, 9);

		UpdateValueProvider<String> move = old -> "MOVED".equals(old) ? old : "MOVED";
		for (int i = 0; i < q.length; i++) {
			a.update(q[i], move);
		}
		RadixKMerStore.BatchBuffers buf = new RadixKMerStore.BatchBuffers(1024);
		for (int i = 0; i < q.length; i++) {
			if (buf.add(q[i])) {
				b.updateBatch(buf, move);
			}
		}
		if (!buf.isEmpty()) {
			b.updateBatch(buf, move);
		}

		int mismatches = 0;
		for (long kmer : kmers) {
			String va = a.getLong(kmer, null);
			String vb = b.getLong(kmer, null);
			if (va == null ? vb != null : !va.equals(vb)) {
				mismatches++;
			}
		}
		if (mismatches != 0) {
			throw new AssertionError("Batched update differs from sequential in " + mismatches + " entries.");
		}
		System.out.printf("correctness OK: batched == sequential over %d entries%n", n);
	}

	private static RadixKMerStore<String> buildStore(long[] kmers) {
		int radixSize = 1 << RADIX_BITS;
		int[] bucketSizes = new int[radixSize];
		for (long kmer : kmers) {
			bucketSizes[RadixKMerStore.radixOf(kmer, RADIX_BITS)]++;
		}
		RadixKMerStore<String> store = new RadixKMerStore<>(K, RADIX_BITS, bucketSizes, 1e-4, 0.01, null, true);
		// putLong resolves values with a lock-free read, so register the whole value pool up front.
		for (int v = 0; v < VALUE_POOL; v++) {
			store.getAddValueIndex(Integer.toString(v));
		}
		Random rnd = new Random(1234);
		for (long kmer : kmers) {
			store.putLong(kmer, Integer.toString(rnd.nextInt(VALUE_POOL)));
		}
		store.optimize();
		return store;
	}

	private static long[] distinctKmers(int n, long seed) {
		Random rnd = new Random(seed);
		long[] out = new long[n];
		for (int i = 0; i < n; i++) {
			out[i] = rnd.nextLong() & KMER_MASK;
		}
		return out;
	}

	private static long[] buildQueries(long[] kmers, int q, long seed) {
		Random rnd = new Random(seed);
		long[] out = new long[q];
		for (int i = 0; i < q; i++) {
			if (rnd.nextInt(10) == 0) {
				out[i] = rnd.nextLong() & KMER_MASK; // likely a miss
			} else {
				out[i] = kmers[rnd.nextInt(kmers.length)]; // a hit
			}
		}
		return out;
	}
}
