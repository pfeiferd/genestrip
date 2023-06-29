package org.metagene.genestrip.trie;

import java.io.Serializable;

import org.metagene.genestrip.util.CGAT;
import org.metagene.genestrip.util.CGATRingBuffer;

public class K31MerSortedArray<V extends Serializable> implements KMerStore<V> {
	private long[][] kmers;	
	private long entries;
	
	public K31MerSortedArray() {
		
	}

	@Override
	public int getLen() {
		return 31;
	}

	@Override
	public long getEntries() {
		return entries;
	}

	@Override
	public void put(CGATRingBuffer buffer, V value, boolean reverse) {
		long kmer = reverse ? CGAT.kmerToLongReverse(buffer) : CGAT.kmerToLongStraight(buffer);
	}

	@Override
	public void put(byte[] nseq, int start, V value, boolean reverse) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public V get(CGATRingBuffer buffer, boolean reverse) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public V get(byte[] nseq, int start, boolean reverse) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void visit(KMerStoreVisitor<V> visitor, boolean reverse) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void optimize() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public boolean isOptimized() {
		// TODO Auto-generated method stub
		return false;
	}

}
