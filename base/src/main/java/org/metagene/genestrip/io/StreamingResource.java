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

import org.metagene.genestrip.util.progressbar.GSProgressUpdate;

import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;

/**
 * A named source of data that can be opened for streaming reads. Opening yields a
 * {@link StreamAccess} that exposes the input stream and tracks how many bytes have been read, so
 * progress can be reported.
 */
public interface StreamingResource {
	/**
	 * A handle on an opened {@link StreamingResource}: it provides the input stream to read from, the
	 * number of bytes read so far and, when known, the total size, and it drives progress updates and
	 * closes the underlying stream.
	 */
	interface StreamAccess extends Closeable, GSProgressUpdate {
		/**
		 * Returns the input stream to read the resource's data from.
		 *
		 * @return the input stream
		 * @throws IOException if the stream cannot be provided
		 */
		public InputStream getInputStream() throws IOException;

		/**
		 * Returns the number of bytes read from the resource so far.
		 *
		 * @return the number of bytes read
		 */
		public long getBytesRead();

		/**
		 * Returns the total size of the resource in bytes, or {@code -1} if unknown.
		 *
		 * @return the total size in bytes, or {@code -1} if unknown
		 * @throws IOException if the size cannot be determined
		 */
		default long getSize() throws IOException {
			return -1;
		}
		
		@Override
		default void close() throws IOException {
			InputStream is = getInputStream();
			if (is != null) {
				is.close();
			}
		}

		@Override
		default long current() {
			return getBytesRead();
		}

		default long max() {
			try {
				return getSize();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
	}
	
	/**
	 * Returns the total size of this resource in bytes, or {@code -1} if unknown.
	 *
	 * @return the total size in bytes, or {@code -1} if unknown
	 * @throws IOException if the size cannot be determined
	 */
	default long getSize() throws IOException {
		return -1;
	}

	/**
	 * Returns the name of this resource.
	 *
	 * @return the resource name
	 */
	public String getName();

	/**
	 * Returns an optional hint about the resource's content type, or {@code null} if none.
	 *
	 * @return the content type hint, or {@code null} if none
	 */
	default String getTypeHint() {
		return null;
	}

	/**
	 * Opens this resource for streaming and returns a handle for reading from it.
	 *
	 * @return a handle for reading from the opened resource
	 * @throws IOException if the resource cannot be opened
	 */
	public StreamAccess openStream() throws IOException;
}
