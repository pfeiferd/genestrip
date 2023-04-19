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
package org.metagene.genestrip.bloom.attic;

import java.io.Serializable;

import org.metagene.genestrip.bloom.TwoLongsCGATBloomFilter;

public class SimpleCGATBloomFilter extends TwoLongsCGATBloomFilter implements Serializable {
	private static final long serialVersionUID = 1L;
	
	private final int[] hashFactors;

	public SimpleCGATBloomFilter(int k, long expectedInsertions, double fpp) {
		super(k, expectedInsertions, fpp);
		this.hashFactors = primeNumbersBruteForce(hashes);
	}

	@Override
	protected int combineLongHashes(long hash1, long hash2, int i) {
		return (int) (hash1 * hashFactors[i] + hash2);
	}

	// 2 is not in it...
	public static int[] primeNumbersBruteForce(int n) {
		int[] res = new int[n];
		int count = 0;
		int num = 3;
		while (count < n) {
			boolean prime = true;// to determine whether the number is prime or not
			int max = (int) Math.sqrt(num);
			for (int i = 2; i <= max; i++) {
				if (num % i == 0) {
					prime = false;
					break;
				}
			}
			if (prime) {
				res[count++] = num;
			}
			num++;
		}
		return res;
	}
}