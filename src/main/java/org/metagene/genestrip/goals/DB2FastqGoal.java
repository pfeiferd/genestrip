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
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.GSFileType;
import org.metagene.genestrip.fastqgen.KMerFastqGenerator;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

public class DB2FastqGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<SmallTaxIdNode>, GSProject> taxNodesGoal;
	private final ObjectGoal<Database, GSProject> storeGoal;
	private final Map<File, SmallTaxIdNode> fileToTaxid;

	@SafeVarargs
	public DB2FastqGoal(GSProject project, ObjectGoal<Set<SmallTaxIdNode>, GSProject> taxNodesGoal,
			ObjectGoal<Database, GSProject> storeGoal, Goal<GSProject>... deps) {
		super(project, GSGoalKey.DB2FASTQ, (List<File>) null, true, Goal.append(deps, taxNodesGoal, storeGoal));
		this.taxNodesGoal = taxNodesGoal;
		this.storeGoal = storeGoal;
		this.fileToTaxid = new HashMap<File, SmallTaxIdNode>();
	}

	@Override
	protected void provideFiles() {
		for (SmallTaxIdNode node : taxNodesGoal.get()) {
			File file = getOutputFile(node.getTaxId());
			addFile(file);
			fileToTaxid.put(file, node);
		}
	}

	protected File getOutputFile(String taxid) {
		return getProject().getOutputFile(getKey().getName(), taxid, null, GSFileType.FASTQ_RES, true);
	}

	@Override
	protected void makeFile(File file) {
		KMerFastqGenerator generator = new KMerFastqGenerator(storeGoal.get().getKmerStore());
		SmallTaxIdNode node = fileToTaxid.get(file);
		try {
			generator.generateFastq(file, node.getTaxId(), getProject().getName() + ":" + node.getName());
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
