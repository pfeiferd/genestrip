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
package org.metagene.genestrip.tax;

import java.util.Collections;
import java.util.Set;

import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * What a {@code taxids.txt} amounts to: the tax ids a database is to hold, and the ones it was told
 * to leave out.
 * <p>
 * The two travel together because they are read together and mean nothing apart. The selected set is
 * what almost everything wants -- it is the database's content, the excluded ones having been
 * subtracted from it already. The excluded set matters to one reader only, the least common ancestor
 * update under
 * {@link org.metagene.genestrip.GSConfigKey.UpdateScope#ALL_BUT_EXCLUDED}, which skips the release's
 * contigs of exactly these taxa. Leaving a branch out of a database and leaving it out of the
 * reckoning are two different things, and this pair is what lets a goal say both at once instead of
 * offering the second through a side door.
 * <p>
 * Both sets are completed down the tree, but not alike: the selected ones stop at the configured
 * rank, the excluded ones do not, since leaving a branch out has to mean the whole branch.
 */
public class TaxNodeSelection {
	private final Set<TaxIdNode> selected;
	private final Set<TaxIdNode> excluded;

	/**
	 * Creates the pair.
	 *
	 * @param selected the tax ids the database is to hold, with the excluded ones already subtracted
	 * @param excluded the tax ids {@code taxids.txt} struck out, with everything below them
	 */
	public TaxNodeSelection(Set<TaxIdNode> selected, Set<TaxIdNode> excluded) {
		this.selected = Collections.unmodifiableSet(selected);
		this.excluded = Collections.unmodifiableSet(excluded);
	}

	/**
	 * Returns the tax ids the database is to hold.
	 *
	 * @return the selected tax ids, with their descendants and without the excluded ones
	 */
	public Set<TaxIdNode> getSelected() {
		return selected;
	}

	/**
	 * Returns the tax ids that were struck out.
	 *
	 * @return the excluded tax ids, with everything below them; empty where nothing was struck out
	 */
	public Set<TaxIdNode> getExcluded() {
		return excluded;
	}

	@Override
	public String toString() {
		return selected.size() + " tax ids, " + excluded.size() + " excluded";
	}
}
