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

import java.io.Serializable;

public class TaxIdInfo implements Serializable, Comparable<TaxIdInfo> {
	private static final long serialVersionUID = 1L;

	public static String ARTIFICIAL_PREFIX = "00";

	public final String taxId;
	protected String name;
	protected short rank;
	protected int position;

	public TaxIdInfo(String taxId, Rank rank) {
		this(taxId, (short) (rank == null ? -1 : rank.ordinal()));
	}

	protected TaxIdInfo(String taxId, short rank) {
		this.taxId = taxId;
		this.rank = rank;		
	}
		
	public int getPosition() {
		return position;
	}
	
	public final Rank getRank() {
		return Rank.byOrdinal(rank);
	}

	public final int getRankOrdinal() {
		return rank;
	}

	public String getName() {
		return name;
	}

	public final String getTaxId() {
		return taxId;
	}

	@Override
	public final int compareTo(TaxIdInfo o) {
		return getPosition() - o.getPosition();
	}

	@Override
	public String toString() {
		return "Node: " + taxId;
	}

	public boolean isArtificialTaxIdInfo() {
		return taxId.charAt(0) == '0' && taxId.charAt(1) == '0';
	}
}
