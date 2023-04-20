 /* “Commons Clause” License Condition v1.0
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

import java.util.Random;

import org.apache.commons.codec.digest.MurmurHash3;

public class CGATMurmurBloomFilter extends TwoLongsCGATBloomFilter {
	private static final long serialVersionUID = 1L;
	private final int[] hashFactors;

	public CGATMurmurBloomFilter(int k, long expectedInsertions, double fpp) {
		super(k, expectedInsertions, fpp);
		Random r = new Random(42);
		hashFactors = new int[hashes];
		for (int i = 0; i < hashFactors.length; i++) {
			hashFactors[i] = r.nextInt();
		}
	}
	
	@Override
	protected int combineLongHashes(long hash1, long hash2, int i) {
		return MurmurHash3.hash32(hash1, hash2, hashFactors[i]);
	}
}