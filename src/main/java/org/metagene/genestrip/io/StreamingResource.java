package org.metagene.genestrip.io;

import java.io.IOException;

import org.metagene.genestrip.io.StreamProvider.ByteCountingInputStreamAccess;

public interface StreamingResource {
	public String getName();
	
	public long getSize() throws IOException;

	public boolean isExists();

	public ByteCountingInputStreamAccess getStreamAccess() throws IOException;
}
