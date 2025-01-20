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
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class StreamingResourceListStream implements StreamingResourceStream {
	private final List<StreamingResource> list;
	
	public StreamingResourceListStream(StreamingResource elem) {
		this(Collections.singletonList(elem));
	}
	
	public StreamingResourceListStream() {
		this(new ArrayList<>());
	}
	
	public StreamingResourceListStream(List<StreamingResource> list) {
		this.list = list;
	}
	
	public List<StreamingResource> getList() {
		return list;
	}
	
	@Override
	public Iterator<StreamingResource> iterator() {
		return list.iterator();
	}
	
	@Override
	public int size() {
		return list.size();
	}
	
	@Override
	public long getTotalByteSize() throws IOException {
		long sum = 0;
		for (StreamingResource r : list) {
			long s = r.getSize();
			if (s < 0) {
				return -1;
			}
			sum += s;
		}
		return sum;
	}

	@Override
	public String toString() {
		return "StreamingResourceListStream:" + (list == null ? "null" : list.toString());
	}
}
