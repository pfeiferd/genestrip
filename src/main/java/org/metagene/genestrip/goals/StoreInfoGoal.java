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
import java.io.PrintStream;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.ResultReporter;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.tax.TaxTree;

public class StoreInfoGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final ObjectGoal<KMerStoreWrapper, GSProject> storeGoal;

	@SafeVarargs
	public StoreInfoGoal(GSProject project, String name, ObjectGoal<TaxTree, GSProject> taxTreeGoal, ObjectGoal<KMerStoreWrapper, GSProject> storeGoal, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.CSV, false),
				Goal.append(deps, taxTreeGoal, storeGoal));
		this.taxTreeGoal = taxTreeGoal;
		this.storeGoal = storeGoal;
	}

	@Override
	protected void makeFile(File file) {
		try {
			KMerStoreWrapper wrapper = storeGoal.get();
			PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
			new ResultReporter(taxTreeGoal.get()).printStoreInfo(wrapper.getStats(), out);
			out.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
