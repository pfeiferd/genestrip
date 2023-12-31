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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
		Map<File, TaxIdNode> res = new HashMap<File, TaxTree.TaxIdNode>();
		File additonalEntryFile = getProject().getAddtionalFile();
		if (additonalEntryFile.exists()) {
			Map<String, List<File>> map;
			map = MultiMatchGoal.readMultiCSV(getProject().getFastaDir(), additonalEntryFile, getLogger());
			TaxTree taxTree = taxTreeGoal.get();
			for (String key : map.keySet()) {
				TaxIdNode node = taxTree.getNodeByTaxId(key);
				if (node != null) {
					for (File file : map.get(key)) {
						res.put(file, node);
					}
				} else {
					if (getLogger().isWarnEnabled()) {
						getLogger().warn("Unknown taxid in additional file (omitting fasta files for it): " + key);
					}
				}
			}
		}
		set(res);
	}
}
