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
package org.metagene.genestrip.goals.refseq;

import java.io.IOException;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;

public class BloomIndexedGoal extends ObjectGoal<MurmurCGATBloomFilter, GSProject> {
	private final FileGoal<GSProject> bloomIndexGoal;

	@SafeVarargs
	public BloomIndexedGoal(GSProject project, String name, FileGoal<GSProject> bloomIndexGoal,
			Goal<GSProject>... dependencies) {
		super(project, name, Goal.append(dependencies, bloomIndexGoal));
		this.bloomIndexGoal = bloomIndexGoal;
	}

	@Override
	public void makeThis() {
		try {
			MurmurCGATBloomFilter filter = MurmurCGATBloomFilter.load(bloomIndexGoal.getFile());
			set(filter);
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}		
	}
	
	void setFilter(MurmurCGATBloomFilter object) {
		set(object);
	}
}
