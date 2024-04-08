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
