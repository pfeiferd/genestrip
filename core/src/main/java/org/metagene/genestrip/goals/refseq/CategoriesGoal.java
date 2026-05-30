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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.HashSet;
import java.util.Set;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.RefSeqCategory;

public class CategoriesGoal<P extends GSProject> extends ObjectGoal<Set<RefSeqCategory>, P> {
	@SafeVarargs
	public CategoriesGoal(P project, Goal<P>... deps) {
		super(project, GSGoalKey.CATEGORIES, deps);
	}

	@Override
	protected void doMakeThis() {
		try {
			set(readFromFile(getProject().getCategoriesFile()));
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected Set<RefSeqCategory> readFromFile(File file) throws IOException {
		Set<RefSeqCategory> res = new HashSet<RefSeqCategory>();

		try (BufferedReader br = new BufferedReader(
				new InputStreamReader(StreamProvider.getInputStreamForFile(file), StandardCharsets.UTF_8))) {
			String line = null;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				// Ignore whole line comments.
				if (!line.startsWith("#")) {
					// Ignore end of line comments.
					int commentIndex = line.indexOf('#');
					if (commentIndex != -1) {
						line = line.substring(0, commentIndex);
					}
					String category = line.trim();
					RefSeqCategory cat = RefSeqCategory.fromDirectoryString(category);
					res.add(cat);
				}
			}
		}

		if (getLogger().isDebugEnabled()) {
			getLogger().debug("Categories for store: " + res);
		}

		return res;
	}
}
