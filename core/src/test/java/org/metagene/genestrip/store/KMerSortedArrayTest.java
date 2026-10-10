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

import java.io.Serializable;

/**
 * Runs the {@link KMerStore} contract of {@link AbstractKMerStoreTest} against
 * {@link KMerSortedArray}, the binary-search store. The store exists so that its lookup path can be
 * measured against {@link RadixKMerStore}'s, and a comparison is only worth making if both answer
 * the same way, which is what this test establishes.
 */
public class KMerSortedArrayTest extends AbstractKMerStoreTest {
	/**
	 * Smaller than the inherited default, and that is the point: the shared test gives every entry a
	 * value of {@code i % } the radix store's value capacity, which is 524,288, whereas this store
	 * indexes its values with a {@code short} and holds {@link KMerSortedArray#MAX_VALUES} of them.
	 * Sixty thousand entries stay below that and still exercise both array layouts and the search.
	 */
	public KMerSortedArrayTest() {
		testSize = 60 * 1000;
		negativeTestSize = testSize;
	}

	@Override
	public <V extends Serializable> KMerStore<V> createKMerStore(Class<V> clazz, int k, long[] kmers) {
		// The same fill-time and post-optimize false-positive probabilities the radix test uses, so
		// neither store is measured against a filter the other does not have.
		return new KMerSortedArray<V>(k, 0.000001, 0.000001, null, false, true, kmers.length);
	}
}
