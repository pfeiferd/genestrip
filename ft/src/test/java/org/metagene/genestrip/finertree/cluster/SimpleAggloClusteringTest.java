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
package org.metagene.genestrip.finertree.cluster;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Test;
import org.metagene.genestrip.finertree.cluster.SimpleAggloClustering.Method;

/**
 * Tests the agglomerative clustering that turns the pairwise similarity of two species into the
 * dendrogram a refinement is built from. The similarities are given by hand and small enough to work
 * the expected merges out on paper, which is what makes the linkage strategies distinguishable:
 * they differ only in the similarity a merged cluster is credited with, and hence in the order the
 * merges happen at all.
 */
public class SimpleAggloClusteringTest {
	/** A similarity given as a full matrix, as the goal computes it from Jaccard indices. */
	private static Similarity matrix(double[][] values) {
		return new Similarity() {
			@Override
			public int values() {
				return values.length;
			}

			@Override
			public double getSimilarity(int i, int j) {
				return values[i][j];
			}
		};
	}

	/** Collects the value indices below the given node, i.e. the members of that cluster. */
	private static Set<Integer> leavesOf(DendrogramNode node) {
		Set<Integer> leaves = new HashSet<>();
		node.visit(new DendrogramNode.Visitor() {
			@Override
			public void preNode(DendrogramNode each) {
				if (each.getChild1() == null && each.getChild2() == null) {
					leaves.add(each.getValueIndex());
				}
			}

			@Override
			public void postNode(DendrogramNode each) {
			}
		});
		return leaves;
	}

	/** Returns every cluster of the tree, the leaves included. */
	private static List<DendrogramNode> allNodes(DendrogramNode root) {
		List<DendrogramNode> nodes = new ArrayList<>();
		root.visit(new DendrogramNode.Visitor() {
			@Override
			public void preNode(DendrogramNode each) {
				nodes.add(each);
			}

			@Override
			public void postNode(DendrogramNode each) {
			}
		});
		return nodes;
	}

	/**
	 * A single species has nothing to be clustered with, and the result must still be a tree - a leaf
	 * carrying that one index.
	 */
	@Test
	public void testOneValueYieldsALeaf() {
		DendrogramNode root = new SimpleAggloClustering(Method.SINGLE_LINKAGE).cluster(matrix(new double[][] { { 1 } }));
		assertNotNull(root);
		assertNull("a single value cannot have been merged with anything", root.getChild1());
		assertEquals(0, root.getValueIndex());
		assertEquals("a lone leaf is the whole tree", 1, root.size());
	}

	/**
	 * Whatever the strategy and whatever the similarities, every species must end up in the tree
	 * exactly once: a refinement that dropped one would lose its k-mers, and one that held it twice
	 * would place them under two nodes.
	 */
	@Test
	public void testEveryValueAppearsExactlyOnce() {
		double[][] values = { { 1.0, 0.9, 0.1, 0.2, 0.15 }, { 0.9, 1.0, 0.12, 0.3, 0.1 },
				{ 0.1, 0.12, 1.0, 0.8, 0.7 }, { 0.2, 0.3, 0.8, 1.0, 0.75 }, { 0.15, 0.1, 0.7, 0.75, 1.0 } };
		for (Method method : Method.values()) {
			DendrogramNode root = new SimpleAggloClustering(method).cluster(matrix(values));
			assertEquals(method + ": every value must be a leaf", 5, leavesOf(root).size());
			// size() counts leaves and inner nodes alike, and merging joins two clusters at a time, so
			// a tree over n values has 2n-1 nodes. Anything else would mean a merge that was not binary.
			assertEquals(method + ": the tree should be binary over all 5 values", 2 * 5 - 1, root.size());
			assertEquals(method.toString(), new HashSet<>(Arrays.asList(0, 1, 2, 3, 4)), leavesOf(root));
		}
	}

	/**
	 * The two most similar species must be merged first, and hence form a cluster of their own,
	 * whichever strategy is used: the strategies differ in what a merged cluster is worth afterwards,
	 * not in which pair looks closest at the start.
	 */
	@Test
	public void testMostSimilarPairIsMergedFirst() {
		// 0 and 3 are far the closest; the rest is deliberately unremarkable.
		double[][] values = { { 1.0, 0.2, 0.1, 0.95 }, { 0.2, 1.0, 0.3, 0.15 }, { 0.1, 0.3, 1.0, 0.2 },
				{ 0.95, 0.15, 0.2, 1.0 } };
		for (Method method : Method.values()) {
			DendrogramNode root = new SimpleAggloClustering(method).cluster(matrix(values));
			boolean found = false;
			for (DendrogramNode node : allNodes(root)) {
				if (leavesOf(node).equals(new HashSet<>(Arrays.asList(0, 3)))) {
					found = true;
					assertEquals(method + ": the cluster must carry the similarity it was merged at", 0.95,
							node.getSimilarity(), 1e-9);
				}
			}
			assertTrue(method + ": the closest pair should form a cluster of its own", found);
		}
	}

	/**
	 * Two groups that are tight within and distant between must come out as two clusters, which is
	 * the case a refinement exists for: the species of one group share their k-mers and are told
	 * apart from the other group.
	 */
	@Test
	public void testTwoSeparateGroupsAreRecovered() {
		// {0,1,2} share almost everything, {3,4} do too, and across the groups there is little.
		double[][] values = { { 1.0, 0.90, 0.88, 0.05, 0.04 }, { 0.90, 1.0, 0.91, 0.03, 0.06 },
				{ 0.88, 0.91, 1.0, 0.05, 0.05 }, { 0.05, 0.03, 0.05, 1.0, 0.93 },
				{ 0.04, 0.06, 0.05, 0.93, 1.0 } };
		for (Method method : Method.values()) {
			DendrogramNode root = new SimpleAggloClustering(method).cluster(matrix(values));
			assertNotNull(method.toString(), root.getChild1());
			Set<Integer> left = leavesOf(root.getChild1());
			Set<Integer> right = leavesOf(root.getChild2());
			Set<Integer> first = new HashSet<>(Arrays.asList(0, 1, 2));
			Set<Integer> second = new HashSet<>(Arrays.asList(3, 4));
			assertTrue(method + ": the two groups should be the two branches below the root, but were "
					+ left + " and " + right,
					(left.equals(first) && right.equals(second)) || (left.equals(second) && right.equals(first)));
		}
	}

	/**
	 * Single linkage merges a cluster with whatever is closest to any of its members, complete
	 * linkage with whatever is closest to all of them. Here that decides where the third species
	 * goes, and the two strategies must disagree - if they did not, the choice would be idle.
	 */
	@Test
	public void testSingleAndCompleteLinkageDisagreeWhereTheyMust() {
		// 0 and 1 merge first. 2 is close to 1 but distant from 0; 3 is moderately close to both.
		// Single linkage credits {0,1} with max(0.10, 0.80) = 0.80 against 2, so 2 joins.
		// Complete linkage credits it with min(0.10, 0.80) = 0.10 against 2 and 0.40 against 3.
		double[][] values = { { 1.0, 0.90, 0.10, 0.40 }, { 0.90, 1.0, 0.80, 0.45 }, { 0.10, 0.80, 1.0, 0.20 },
				{ 0.40, 0.45, 0.20, 1.0 } };
		Set<Integer> zeroOneTwo = new HashSet<>(Arrays.asList(0, 1, 2));
		Set<Integer> zeroOneThree = new HashSet<>(Arrays.asList(0, 1, 3));

		boolean singleJoinsTwo = false;
		for (DendrogramNode node : allNodes(new SimpleAggloClustering(Method.SINGLE_LINKAGE).cluster(matrix(values)))) {
			singleJoinsTwo |= leavesOf(node).equals(zeroOneTwo);
		}
		boolean completeJoinsThree = false;
		for (DendrogramNode node : allNodes(new SimpleAggloClustering(Method.COMPLETE_LINKAGE).cluster(matrix(values)))) {
			completeJoinsThree |= leavesOf(node).equals(zeroOneThree);
		}
		assertTrue("single linkage should let 2 join {0,1}, as it is close to 1", singleJoinsTwo);
		assertTrue("complete linkage should let 3 join {0,1}, as 2 is distant from 0", completeJoinsThree);
	}

	/**
	 * Merging proceeds from the most similar pair to the least, so a cluster can never be more
	 * similar within than the one it was merged into. A dendrogram that broke this could not be cut
	 * at a similarity threshold, which is how the refinement reads it.
	 */
	@Test
	public void testSimilarityDoesNotIncreaseTowardsTheRoot() {
		double[][] values = { { 1.0, 0.90, 0.30, 0.20, 0.10 }, { 0.90, 1.0, 0.35, 0.25, 0.15 },
				{ 0.30, 0.35, 1.0, 0.85, 0.40 }, { 0.20, 0.25, 0.85, 1.0, 0.45 },
				{ 0.10, 0.15, 0.40, 0.45, 1.0 } };
		for (Method method : Method.values()) {
			DendrogramNode root = new SimpleAggloClustering(method).cluster(matrix(values));
			assertMonotone(root, method);
		}
	}

	private static void assertMonotone(DendrogramNode node, Method method) {
		for (DendrogramNode child : new DendrogramNode[] { node.getChild1(), node.getChild2() }) {
			if (child != null) {
				assertTrue(method + ": a child cluster (" + child.getSimilarity() + ") should not be less"
						+ " similar than its parent (" + node.getSimilarity() + ")",
						child.getSimilarity() >= node.getSimilarity() - 1e-9);
				assertMonotone(child, method);
			}
		}
	}
}
