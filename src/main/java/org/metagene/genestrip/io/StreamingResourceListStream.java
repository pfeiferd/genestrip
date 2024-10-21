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
}
