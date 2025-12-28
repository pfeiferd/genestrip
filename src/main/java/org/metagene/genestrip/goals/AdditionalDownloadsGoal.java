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
import java.io.InputStreamReader;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.Goal;

public class AdditionalDownloadsGoal<P extends GSProject> extends GSFileDownloadGoal<P> {
	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#')
			.setDelimiter(' ').setRecordSeparator('\n').build();

	public static final String DOWNLOADS_NAME = "downloads.txt";

	private List<File> files;
	private final Map<File, URL> fileToURL;
	private final Map<File, String> fileToMD5;

	@SafeVarargs
	public AdditionalDownloadsGoal(P project, Goal<P>... deps) {
		super(project, GSGoalKey.ADD_DOWNLOADS, deps);
		fileToURL = new HashMap<File, URL>();
		fileToMD5 = new HashMap<File, String>();
	}

	@Override
	public boolean isAllowTransitiveClean() {
		return false;
	}

	@Override
	public List<File> getFiles() {
		if (files == null) {
			files = new ArrayList<File>();
			File additonalEntryFile = getProject().getAdditionalFile();
			if (additonalEntryFile.exists()) {
				List<String> forGlobalDownload = new ArrayList<String>();
				try (CSVParser parser = FORMAT
						.parse(new InputStreamReader(StreamProvider.getInputStreamForFile(additonalEntryFile)))) {
					for (CSVRecord record : parser) {
						String fileName = record.get(1);
						if (record.size() >= 3) {
							String urlString = record.get(2);
							String md5 = record.size() >= 4 ? record.get(3) : null;
							File file = new File(getProject().getFastaDir(), fileName);
							URL url = new URL(urlString);
							fileToURL.put(file, url);
							if (md5 != null) {
								fileToMD5.put(file, md5);
							}
							if (getLogger().isDebugEnabled()) {
								getLogger().debug("Additional file: " + file);
							}
							files.add(file);
						} else {
							if (!new File(getProject().getFastaDir(), fileName).exists()) {
								forGlobalDownload.add(fileName);
							}
						}
					}
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
				File downloadsFile = new File(getProject().getCommon().getFastaDir(), DOWNLOADS_NAME);
				if (downloadsFile.exists()) {
					try (CSVParser parser = FORMAT
							.parse(new InputStreamReader(StreamProvider.getInputStreamForFile(downloadsFile)))) {
						for (CSVRecord record : parser) {
							if (record.size() >= 2) {
								String fileName = record.get(0);
								if (forGlobalDownload.contains(fileName)) {
									String urlString = record.get(1);
									String md5 = record.size() >= 3 ? record.get(2) : null;
									if (urlString != null) {
										File file = new File(getProject().getCommon().getFastaDir(), fileName);
										URL url = new URL(urlString);

										fileToURL.put(file, url);
										if (md5 != null) {
											fileToMD5.put(file, md5);
										}
										files.add(file);

										if (getLogger().isDebugEnabled()) {
											getLogger().debug("Additional download file: " + file);
										}
									}
								}
							}
						}
					} catch (IOException e) {
						throw new RuntimeException(e);
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

		if (getLogger().isDebugEnabled()) {
			getLogger().debug("Additional download for " + url.toExternalForm());
			getLogger().debug("Saving file " + file.toString());
		}
		try (ReadableByteChannel readableByteChannel = Channels.newChannel(url.openStream());
				FileOutputStream out = new FileOutputStream(file)) {
			out.getChannel().transferFrom(readableByteChannel, 0, Long.MAX_VALUE);
		}
	}

	@Override
	protected String getMD5CheckSum(File file) {
		return fileToMD5.get(file);
	}
}
