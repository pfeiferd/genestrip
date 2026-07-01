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

import java.io.IOException;

/**
 * An iterable sequence of {@link StreamingResource}s, e.g. the set of input resources to be
 * processed one after another.
 */
public interface StreamingResourceStream extends Iterable<StreamingResource> {
	/**
	 * Returns the combined size in bytes of all resources in the stream, or {@code -1} if unknown.
	 *
	 * @return the combined size in bytes, or {@code -1} if unknown
	 * @throws java.io.IOException if determining the size requires I/O that fails
	 */
	default long getTotalByteSize() throws IOException {
		return -1;
	}

	/**
	 * Returns the number of resources in the stream, or {@code -1} if unknown.
	 *
	 * @return the number of resources, or {@code -1} if unknown
	 */
	// Number of elements in the stream, maybe -1 if known.
	default int size() {
		return -1;
	}
	
	/**
	 * Unchecked wrapper for an {@link IOException} raised while iterating a
	 * {@link StreamingResourceStream}, allowing it to propagate through the {@link Iterable} API.
	 */
	public static class IteratorRuntimeIOException extends RuntimeException {
		private static final long serialVersionUID = 1L;

		/**
		 * Wraps the given {@link IOException}.
		 *
		 * @param e the I/O exception to wrap
		 */
		public IteratorRuntimeIOException(IOException e) {
			super(e);
		}
	}
}
