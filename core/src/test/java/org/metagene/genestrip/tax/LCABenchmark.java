package org.metagene.genestrip.tax;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * Stand-alone benchmark + correctness check for the O(d1+d2) {@link TaxTree#getLowestCommonAncestor}
 * (optimization 5b) against the previous O(d1*d2) nested-scan reference, on the real NCBI taxonomy.
 *
 * Run (no network): java -cp <cp> org.metagene.genestrip.tax.LCABenchmark [taxDir] [pairs]
 * taxDir defaults to ./data/common (must contain nodes.dmp and names.dmp).
 */
public class LCABenchmark {
	public static void main(String[] args) {
		File taxDir = new File(args.length > 0 ? args[0] : "./data/common");
		int pairs = args.length > 1 ? Integer.parseInt(args[1]) : 5_000_000;

		System.out.println("Loading taxonomy from " + taxDir + " ...");
		TaxTree tree = new TaxTree(taxDir, false);
		List<TaxIdNode> nodes = collectNodes(tree);
		int maxDepth = 0;
		long sumDepth = 0;
		for (TaxIdNode n : nodes) {
			int d = depthOf(n);
			sumDepth += d;
			if (d > maxDepth) {
				maxDepth = d;
			}
		}
		System.out.printf("nodes=%d maxDepth=%d avgDepth=%.1f%n", nodes.size(), maxDepth, (double) sumDepth / nodes.size());

		// Correctness: new impl must match the reference on a large random sample (incl. equal nodes,
		// ancestor pairs and fully divergent pairs).
		Random rnd = new Random(1);
		int n = nodes.size();
		for (int i = 0; i < 2_000_000; i++) {
			TaxIdNode a = nodes.get(rnd.nextInt(n));
			TaxIdNode b = rnd.nextInt(8) == 0 ? a : nodes.get(rnd.nextInt(n));
			if (referenceLCA(a, b) != tree.getLowestCommonAncestor(a, b)) {
				throw new AssertionError("LCA mismatch for " + a.getTaxId() + " / " + b.getTaxId());
			}
		}
		System.out.println("correctness OK: new LCA matches reference over 2,000,000 pairs");

		// Two workloads: fully random (divergent, the nested loop's worst case) and "close" pairs -
		// node2 is a descendant of one of node1's near ancestors, i.e. the ancestor/nearby case that
		// dominates the real update (a k-mer stored at genus level found in a species of that genus).
		runWorkload("random  ", tree, nodes, randomPairs(rnd, n, pairs, false));
		runWorkload("close   ", tree, nodes, randomPairs(rnd, n, pairs, true));
	}

	private static int[][] randomPairs(Random rnd, int n, int pairs, boolean close) {
		int[] ia = new int[pairs];
		int[] ib = new int[pairs];
		for (int i = 0; i < pairs; i++) {
			ia[i] = rnd.nextInt(n);
			ib[i] = close ? ia[i] : rnd.nextInt(n); // "close" resolves node2 near node1 below
		}
		return new int[][] { ia, ib };
	}

	private static void runWorkload(String label, TaxTree tree, List<TaxIdNode> nodes, int[][] idx) {
		int[] ia = idx[0];
		int[] ib = idx[1];
		boolean close = ia == ib || java.util.Arrays.equals(ia, ib);
		// For "close", derive node2 = walk node1 up a few steps then reuse node1's neighbourhood.
		TaxIdNode[] a = new TaxIdNode[ia.length];
		TaxIdNode[] b = new TaxIdNode[ia.length];
		Random rnd = new Random(99);
		for (int i = 0; i < ia.length; i++) {
			a[i] = nodes.get(ia[i]);
			if (close) {
				TaxIdNode up = a[i];
				int steps = rnd.nextInt(4);
				for (int s = 0; s < steps && up.getParent() != null; s++) {
					up = up.getParent();
				}
				b[i] = up; // an ancestor of node1 (0-3 levels up): the common ancestor/close case
			} else {
				b[i] = nodes.get(ib[i]);
			}
		}
		timeRef(a, b);
		timeNw(tree, a, b);
		long refNs = 0, newNs = 0;
		int reps = 3;
		for (int r = 0; r < reps; r++) {
			refNs += timeRef(a, b);
			newNs += timeNw(tree, a, b);
		}
		double refPer = (double) refNs / ((long) reps * a.length);
		double newPer = (double) newNs / ((long) reps * a.length);
		System.out.printf("%n[%s] reference O(d1*d2): %6.1f ns/LCA   new O(d1+d2): %6.1f ns/LCA   speedup %.2fx%n",
				label, refPer, newPer, refPer / newPer);
	}

	private static long timeRef(TaxIdNode[] a, TaxIdNode[] b) {
		long acc = 0;
		long t0 = System.nanoTime();
		for (int i = 0; i < a.length; i++) {
			if (referenceLCA(a[i], b[i]) != null) {
				acc++;
			}
		}
		long dt = System.nanoTime() - t0;
		if (acc == Long.MIN_VALUE) {
			System.out.println(acc);
		}
		return dt;
	}

	private static long timeNw(TaxTree tree, TaxIdNode[] a, TaxIdNode[] b) {
		long acc = 0;
		long t0 = System.nanoTime();
		for (int i = 0; i < a.length; i++) {
			if (tree.getLowestCommonAncestor(a[i], b[i]) != null) {
				acc++;
			}
		}
		long dt = System.nanoTime() - t0;
		if (acc == Long.MIN_VALUE) {
			System.out.println(acc);
		}
		return dt;
	}

	// The previous implementation, kept here as the reference oracle.
	private static TaxIdNode referenceLCA(TaxIdNode node1, TaxIdNode node2) {
		if (node1 == node2) {
			return node1;
		}
		for (TaxIdNode res = node1; res != null; res = res.getParent()) {
			for (TaxIdNode ancestor2 = node2; ancestor2 != null; ancestor2 = ancestor2.getParent()) {
				if (res == ancestor2) {
					return res;
				}
			}
		}
		return null;
	}

	private static int depthOf(TaxIdNode node) {
		int d = 0;
		for (TaxIdNode n = node.getParent(); n != null; n = n.getParent()) {
			d++;
		}
		return d;
	}

	private static List<TaxIdNode> collectNodes(TaxTree tree) {
		List<TaxIdNode> out = new ArrayList<>();
		collect(tree.getRoot(), out);
		return out;
	}

	private static void collect(TaxIdNode node, List<TaxIdNode> out) {
		if (node == null) {
			return;
		}
		out.add(node);
		List<TaxIdNode> subs = node.getSubNodes();
		if (subs != null) {
			for (int i = 0; i < subs.size(); i++) {
				collect(subs.get(i), out);
			}
		}
	}
}
