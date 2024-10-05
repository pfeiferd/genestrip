package org.metagene.genestrip.io;

import java.io.IOException;

public interface StreamingResourceStream extends Iterable<StreamingResource> {
	default long getTotalByteSize() throws IOException {
		return -1;
	}
	
	// Number of elements in the stream, maybe -1 if known.
	default int size() {
		return -1;
	}
	
	public static class IteratorRuntimeIOException extends RuntimeException {
		private static final long serialVersionUID = 1L;

		public IteratorRuntimeIOException(IOException e) {
			super(e);
		}
	}
}
