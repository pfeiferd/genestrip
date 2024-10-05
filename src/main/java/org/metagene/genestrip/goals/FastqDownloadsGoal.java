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
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.io.StreamingResourceListStream;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.io.StreamingURLResource;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;

public class FastqDownloadsGoal extends GSFileDownloadGoal {
	private List<File> files;
	private final Map<File, URL> fileToURL;
	private final ObjectGoal<Map<String, StreamingResourceStream>, GSProject> mapGoal;
	private final ObjectGoal<Map<String, StreamingResourceStream>, GSProject> transformGoal;

	@SafeVarargs
	public FastqDownloadsGoal(GSProject project, ObjectGoal<Map<String, StreamingResourceStream>, GSProject> mapGoal,
			ObjectGoal<Map<String, StreamingResourceStream>, GSProject> transformGoal, Goal<GSProject>... deps) {
		super(project, GSGoalKey.FASTQ_DOWNLOAD, Goal.append(deps, mapGoal, transformGoal));
		fileToURL = new HashMap<File, URL>();
		this.mapGoal = mapGoal;
		this.transformGoal = transformGoal;
	}

	@Override
	public boolean isAllowTransitiveClean() {
		return false;
	}

	@Override
	public List<File> getFiles() {
		if (files == null) {
			files = new ArrayList<File>();
			for (String key : mapGoal.get().keySet()) {
				int index = 0;
				StreamingResourceStream list = mapGoal.get().get(key);
				if (list instanceof StreamingResourceListStream) {
					StreamingResourceStream ll = (StreamingResourceListStream) list;
					for (StreamingResource resource : ll) {
						if (resource instanceof StreamingURLResource) {
							StreamingResource fr = ((StreamingResourceListStream) transformGoal.get().get(key))
									.getList().get(index);
							if (fr instanceof StreamingFileResource) {
								File file = ((StreamingFileResource) fr).getFile();
								files.add(file);
								fileToURL.put(file, ((StreamingURLResource) resource).getURL());
							}
						}
						index++;
					}
				}
			}
		}
		return files;
	}

	@Override
	protected String getFTPDir(File file) {
		throw new IllegalStateException("Should not be called from here...");
	}

	@Override
	protected boolean isAdditionalFile(File file) {
		return true;
	}

	@Override
	public void additionalDownload(File file) throws IOException {
		URL url = fileToURL.get(file);

		if (getLogger().isInfoEnabled()) {
			getLogger().info("Fastq download for " + url.toExternalForm());
			getLogger().info("Saving file " + file.toString());
		}
		try (ReadableByteChannel readableByteChannel = Channels.newChannel(url.openStream());
				FileOutputStream out = new FileOutputStream(file)) {
			out.getChannel().transferFrom(readableByteChannel, 0, Long.MAX_VALUE);
		}
	}
}
