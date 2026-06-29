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

import org.metagene.genestrip.bloom.KMerProbFilter;

/**
 * Optional capability interface for {@link KMerStore} implementations that support a
 * probabilistic pre-filter for lookups.
 * <p>
 * These operations are not part of the core {@link KMerStore} contract because an exchangeable
 * store implementation may not offer them. Callers should therefore feature-check via
 * {@code instanceof TunableKMerStore} before using them.
 */
public interface TunableKMerStore<V extends Serializable> extends KMerStore<V> {
	public KMerProbFilter getFilter();

	public void setFilter(KMerProbFilter filter);

	public void setUseFilter(boolean useFilter);

	public boolean isUseFilter();
}
