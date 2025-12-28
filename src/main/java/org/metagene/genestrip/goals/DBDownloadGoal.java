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
import java.io.InputStream;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.Collections;
import java.util.List;

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.ByteCountingInputStream;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.io.StreamingURLResource;
import org.metagene.genestrip.make.Goal;

public class DBDownloadGoal<P extends GSProject> extends GSFileDownloadGoal<P> {
	private final StreamingResource dbResource;
	private final String md5;
	private List<File> dbFile;
	private boolean dumped;

	@SafeVarargs
	public DBDownloadGoal(P project, StreamingResource dbResource, String md5, Goal<P>... deps) {
		super(project, GSGoalKey.DB_DOWNLOAD, deps);
		this.dbResource = dbResource;
		this.md5 = md5;
	}

	@SafeVarargs
	public DBDownloadGoal(P project, URL url, String md5, long logCycle, Goal<P>... deps) {
		super(project, GSGoalKey.DB_DOWNLOAD, deps);
		this.dbResource = new StreamingURLResource(url.getFile(), url, true) {
			@Override
			protected ByteCountingInputStream createByteCountingInputStream(InputStream is) throws IOException {
				return new ByteCountingInputStream(is) {
					@Override
					public int read() throws IOException {
						if (dumped) {
							throw new DBDownloadInterruptedException();
						}
						if (logCycle > 0 && bRead % logCycle == 0) {
							log(bRead);
						}
						return super.read();
					}

					@Override
					public int read(byte[] b) throws IOException {
						if (dumped) {
							throw new DBDownloadInterruptedException();
						}
						if (logCycle > 0 && bRead % logCycle == 0) {
							log(bRead);
						}
						return super.read(b);
					}

					@Override
					public int read(byte[] b, int off, int len) throws IOException {
						if (dumped) {
							throw new DBDownloadInterruptedException();
						}
						if (logCycle > 0 && bRead % logCycle == 0) {
							log(bRead);
						}
						return super.read(b, off, len);
					}
				};
			}
		};
		this.md5 = md5;
	}

	public boolean isDumped() {
		return dumped;
	}

	@Override
	public void dump() {
		super.dump();
		dumped = true;
	}

	protected void log(long bytesCovered) {
	}

	@Override
	public boolean isAllowTransitiveClean() {
		return false;
	}

	@Override
	protected String getMD5CheckSum(File file) {
		return md5;
	}

	@Override
	public List<File> getFiles() {
		if (dbFile == null) {
			dbFile = Collections.singletonList(getProject().getDBFile());
		}
		return dbFile;
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
		if (getLogger().isDebugEnabled()) {
			getLogger().debug("DB download for " + dbResource);
			getLogger().debug("Saving DB file " + file.toString());
		}
		try (StreamingResource.StreamAccess lbyteCountAccess = dbResource.openStream()) {
			try (ReadableByteChannel readableByteChannel = Channels.newChannel(lbyteCountAccess.getInputStream());
					FileOutputStream out = new FileOutputStream(file)) {
				out.getChannel().transferFrom(readableByteChannel, 0, Long.MAX_VALUE);
			}
		}
		if (Thread.interrupted()) {
			dumped = true;
		}
		if (getLogger().isWarnEnabled()) {
			if (dumped) {
				getLogger().warn("Interruption when saving DB file " + file.toString() + ". Saving canceled.");
			} else {
				getLogger().debug("Saved DB file " + file.toString());
			}
		}
	}

	public static class DBDownloadInterruptedException extends RuntimeException {
		private static final long serialVersionUID = 1L;

		public DBDownloadInterruptedException() {
			super();
		}
	}
}
