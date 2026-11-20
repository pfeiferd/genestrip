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

import java.io.File;
import java.io.IOException;

import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * The file format a {@link TaxTree} is read from. An implementation parses the taxonomy files and
 * builds the nodes through the {@link TaxTree.Builder} it is handed; the tree itself knows nothing
 * of any format. {@link NCBITaxTreeFormat} reads the NCBI taxonomy dump.
 */
public interface TaxTreeFormat {
	/**
	 * Reads the taxonomy files in the given directory and creates their nodes, links and names via
	 * the builder.
	 *
	 * @param path        the directory containing the taxonomy files
	 * @param builder     creates and links the nodes of the tree being read
	 * @param progressBar whether to show a progress bar while reading the files
	 * @return the root node
	 * @throws IOException if reading the files fails
	 */
	public TaxIdNode read(File path, TaxTree.Builder builder, boolean progressBar) throws IOException;
}
