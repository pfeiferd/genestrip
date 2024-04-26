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
package org.metagene.genestrip.io;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import org.metagene.genestrip.io.StreamProvider.ByteCountingInputStreamAccess;

public class StreamingFileResource implements StreamingResource {
	private final File file;

	public StreamingFileResource(File file) {
		this.file = file;
	}

	@Override
	public long getSize() throws IOException {
		return Files.size(file.toPath());
	}

	public ByteCountingInputStreamAccess getStreamAccess() throws IOException {
		return StreamProvider.getByteCountingInputStreamForFile(file, false);
	}

	@Override
	public boolean isExists() {
		return file.exists();
	}
	
	@Override
	public String getName() {
		return file.getName();
	}
	
	public File getFile() {
		return file;
	}
	
	@Override
	public String toString() {
		return "Streaming File: " + file.toString();
	}
}
