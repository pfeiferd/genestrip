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
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.zip.GZIPInputStream;

/**
 * A {@link StreamingResource} backed by a {@link URL}. The content is decompressed with GZIP unless
 * suppressed, and the number of bytes read is tracked via a {@link ByteCountingInputStream}.
 */
public class StreamingURLResource implements StreamingResource {
	private final String name;
	private final URL url;
	private final boolean noGZ;

	/**
	 * Creates a resource for the given URL, named after its file part, with GZIP decompression.
	 *
	 * @param url the URL backing the resource
	 */
	public StreamingURLResource(URL url) {
		this(url.getFile(), url, false);
	}

	/**
	 * Creates a resource for the given URL with the given name. If {@code noGZ} is {@code true}, the
	 * content is read without GZIP decompression.
	 *
	 * @param name the resource name
	 * @param url the URL backing the resource
	 * @param noGZ whether to read the content without GZIP decompression
	 */
	public StreamingURLResource(String name, URL url, boolean noGZ) {
		this.name = name;
		this.url = url;
		this.noGZ = noGZ;
	}

	/**
	 * Returns whether the content is read without GZIP decompression.
	 *
	 * @return {@code true} if GZIP decompression is suppressed
	 */
	public boolean isNoGZ() {
		return noGZ;
	}

	public StreamingResource.StreamAccess openStream() throws IOException {
		URLConnection connection = url.openConnection();
		connection.connect();
		ByteCountingInputStream is = createByteCountingInputStream(connection.getInputStream());
		InputStream readFromIs = noGZ ? is : new GZIPInputStream(is, StreamProvider.getBufferSize());

		return new StreamAccess() {
			@Override
			public InputStream getInputStream() {
				return readFromIs;
			}

			@Override
			public long getBytesRead() {
				return is.getBytesRead();
			}

			@Override
			public long getSize() throws IOException {
				return connection.getContentLengthLong();
			}
		};
	}

	/**
	 * Wraps the given stream in a {@link ByteCountingInputStream}; may be overridden to customize the
	 * byte-counting stream.
	 *
	 * @param is the stream to wrap
	 * @return the byte-counting stream wrapping {@code is}
	 * @throws IOException if wrapping the stream fails
	 */
	protected ByteCountingInputStream createByteCountingInputStream(InputStream is) throws IOException {
		return new ByteCountingInputStream(is);
	}

	@Override
	public String getName() {
		return name;
	}

	/**
	 * Returns the URL backing this resource.
	 *
	 * @return the backing URL
	 */
	public URL getURL() {
		return url;
	}

	@Override
	public String toString() {
		return "Streaming URL: " + url.toString();
	}
}
