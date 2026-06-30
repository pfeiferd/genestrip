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

import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.BufferedLineReader;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.util.ByteArrayUtil;

public class CheckSumMapGoal<P extends GSProject> extends ObjectGoal<Map<String, String>, P> {
	private final RefSeqCatalogDownloadGoal catalogGoal;

	@SafeVarargs
	public CheckSumMapGoal(P project, RefSeqCatalogDownloadGoal catalogGoal,
			Goal<P>... deps) {
		super(project, GSGoalKey.CHECKSUMMAP, Goal.append(deps, catalogGoal));
		this.catalogGoal = catalogGoal;
	}

	@Override
	protected void doMakeThis() {
		try (BufferedLineReader lineReader = new BufferedLineReader(
				new FileInputStream(catalogGoal.getInstalledFilesFile()))) {
			byte[] target = new byte[2048];
			Map<String, String> map = new HashMap<String, String>();
			int len;
			while ((len = lineReader.nextLine(target)) > 0) {
				// Trim trailing line terminator(s) ('\n', or '\r\n' for CRLF / a final line without
				// one). Cap at target.length since nextLine returns length+1 when a line overflows.
				int end = Math.min(len, target.length);
				while (end > 0 && (target[end - 1] == '\n' || target[end - 1] == '\r')) {
					end--;
				}
				int tab = ByteArrayUtil.indexOf(target, 0, end, '\t');
				if (tab < 0) {
					continue; // skip malformed lines without a tab
				}
				String md5 = new String(target, 0, tab);
				String file = new String(target, tab + 1, end - tab - 1);
				map.put(file, md5);
			}
			set(map);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	@Override
	protected String logSetObject(Map<String, String> object) {
		return "Map of size " + object.size();
	}
}
