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
package org.metagene.genestrip.goals;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class AdditionalFastasGoal extends ObjectGoal<Map<File, TaxIdNode>, GSProject> {
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;

	@SafeVarargs
	public AdditionalFastasGoal(GSProject project, String name, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			Goal<GSProject>... dependencies) {
		super(project, name, append(dependencies, taxTreeGoal));
		this.taxTreeGoal = taxTreeGoal;
	}

	@Override
	public void makeThis() {
		try {
			Map<File, TaxIdNode> res = new HashMap<File, TaxTree.TaxIdNode>();

			File additonalEntryFile = getProject().getAddtionalFile();
			if (additonalEntryFile.exists()) {
				TaxTree taxTree = taxTreeGoal.get();

				Iterable<CSVRecord> records = MultiMatchGoal.readCSVFile(additonalEntryFile);
				for (CSVRecord record : records) {
					String taxid = record.get(0);
					String fileName = record.get(1);

					TaxIdNode node = taxTree.getNodeByTaxId(taxid);
					if (node != null) {
						File file = new File(fileName);
						if (!file.exists()) {
							file = new File(getProject().getFastaDir(), fileName);
						}
						if (file.exists()) {
							if (getLogger().isInfoEnabled()) {
								getLogger()
										.info("Adding additional fasta file " + file + " for taxid " + node.getTaxId());
							}
							res.put(file, node);
						} else {
							if (getLogger().isWarnEnabled()) {
								getLogger().warn("Ignoring missing additional fasta file " + file);
							}
						}
					}
				}
			}
			set(res);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
