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

public class StreamingURLResource implements StreamingResource {
	private final String name;
	private final URL url;
	private final boolean noGZ;

	public StreamingURLResource(URL url) {
		this(url.getFile(), url, false);
	}

	public StreamingURLResource(String name, URL url, boolean noGZ) {
		this.name = name;
		this.url = url;
		this.noGZ = noGZ;
	}

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

	protected ByteCountingInputStream createByteCountingInputStream(InputStream is) throws IOException {
		return new ByteCountingInputStream(is);
	}

	@Override
	public String getName() {
		return name;
	}

	public URL getURL() {
		return url;
	}

	@Override
	public String toString() {
		return "Streaming URL: " + url.toString();
	}
}
