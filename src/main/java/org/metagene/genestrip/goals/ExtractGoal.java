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

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.GSFileType;
import org.metagene.genestrip.fastq.AbstractLoggingFastqStreamer;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.util.ByteArrayUtil;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Map;

public class ExtractGoal<P extends GSProject> extends Goal<P> {
	private final ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal;
	private final ExecutionContext bundle;

	@SafeVarargs
	public ExtractGoal(P project, ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal,
					   ExecutionContext bundle, Goal<P>... deps) {
		super(project, GSGoalKey.EXTRACT, Goal.append(deps, fastqMapGoal));
		this.fastqMapGoal = fastqMapGoal;
		this.bundle = bundle;
	}

	@Override
	public boolean isMade() {
		return false;
	}

	@Override
	protected void doMakeThis() {
		AbstractLoggingFastqStreamer matcher = null;
		try {
			Map<String, StreamingResourceStream> map = fastqMapGoal.get();
			String filter = stringConfigValue(GSConfigKey.EXTRACT_KEY);
			if (filter == null || filter.isEmpty()) {
				if (getLogger().isDebugEnabled()) {
					getLogger().error("Missing extract key in config entry " + GSConfigKey.EXTRACT_KEY + ".");
				}
				return;
			}

			PrintStream[] psRef = new PrintStream[1];
			for (String key : map.keySet()) {
				StreamingResourceStream fastqs = map.get(key);
				if (matcher == null) {
					matcher = new AbstractLoggingFastqStreamer(intConfigValue(GSConfigKey.KMER_SIZE), intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES),
							intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE), bundle, true) {
						@Override
						protected void nextEntry(ReadEntry readStruct, int threadIndex) throws IOException {
							if (ByteArrayUtil.startsWith(readStruct.readDescriptor, 1, filter)) {
								readStruct.print(psRef[0]);
							}
						}

						@Override
						protected boolean isProgressBar() {
							return booleanConfigValue(GSConfigKey.PROGRESS_BAR);
						}

						@Override
						protected String getProgressBarTaskName() {
							return getKey().getName();
						}
					};
				}
				if (booleanConfigValue(GSConfigKey.WRITE_FILTERED_FASTQ)) {
					File filteredFile = getProject().getOutputFile(getKey().getName(), key, null, GSFileType.FASTQ_RES,
							false);
					try (PrintStream ps = new PrintStream(filteredFile)) {
						psRef[0] = ps;
						matcher.processFastqStreams(fastqs);
					}
				}
				else {
					psRef[0] = System.out;
					matcher.processFastqStreams(fastqs);
				}
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			if (matcher != null) {
				matcher.dump();
			}
		}
	}
}
