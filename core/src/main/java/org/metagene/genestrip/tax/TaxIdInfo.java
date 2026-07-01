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

/**
 * Base class for a node in the taxonomy tree. It carries a tax id, a name, a rank
 * (stored as the rank's enum ordinal) and a {@code position} that reflects the node's
 * order within the tree. Nodes are ordered by that position.
 */
public abstract class TaxIdInfo implements Serializable, Comparable<TaxIdInfo> {
	private static final long serialVersionUID = 1L;

	/** The tax id prefix that marks an artificial (non-NCBI) node. */
	public static String ARTIFICIAL_PREFIX = "00";

	/** The tax id of this node. */
	public final String taxId;
	/** The name of this node. */
	protected String name;
	/** The rank of this node, stored as the rank enum's ordinal. */
	protected short rank;
	/** The position of this node within the taxonomy tree. */
	protected int position;

	/**
	 * Creates a node with the given tax id and rank.
	 *
	 * @param taxId the tax id
	 * @param rank  the rank, or {@code null} for no rank
	 */
	public TaxIdInfo(String taxId, Rank rank) {
		this(taxId, (short) (rank == null ? -1 : rank.ordinal()));
	}

	/**
	 * Creates a node with the given tax id and rank ordinal.
	 *
	 * @param taxId the tax id
	 * @param rank  the rank's enum ordinal, or -1 for no rank
	 */
	protected TaxIdInfo(String taxId, short rank) {
		this.taxId = taxId;
		this.rank = rank;		
	}
		
	/**
	 * Returns the position of this node within the taxonomy tree.
	 *
	 * @return the position within the tree
	 */
	public int getPosition() {
		return position;
	}

	/**
	 * Returns the rank of this node.
	 *
	 * @return the rank
	 */
	public final Rank getRank() {
		return Rank.byOrdinal(rank);
	}

	/**
	 * Returns the rank of this node as its enum ordinal.
	 *
	 * @return the rank ordinal
	 */
	public final int getRankOrdinal() {
		return rank;
	}

	/**
	 * Returns the name of this node.
	 *
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * Returns the tax id of this node.
	 *
	 * @return the tax id
	 */
	public final String getTaxId() {
		return taxId;
	}

	/**
	 * Orders nodes by their {@code position} within the taxonomy tree.
	 */
	@Override
	public final int compareTo(TaxIdInfo o) {
		return getPosition() - o.getPosition();
	}

	@Override
	public String toString() {
		return "Node: " + taxId;
	}

	/**
	 * Whether this is an artificial node, i.e. one whose tax id starts with the
	 * {@code "00"} prefix rather than being a real NCBI tax id.
	 *
	 * @return {@code true} if this is an artificial node
	 */
	public boolean isArtificialTaxIdInfo() {
		return taxId.charAt(0) == '0' && taxId.charAt(1) == '0';
	}

	/**
	 * Returns the parent node of this node.
	 *
	 * @return the parent node, or {@code null} if this is the root
	 */
	public abstract TaxIdInfo getParent();
}
