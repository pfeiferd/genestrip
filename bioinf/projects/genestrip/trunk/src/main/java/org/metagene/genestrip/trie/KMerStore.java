package org.metagene.genestrip.trie;

import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import org.metagene.genestrip.util.CGATRingBuffer;
import org.metagene.genestrip.util.StreamProvider;

public interface KMerStore<V extends Serializable> extends Serializable {
	public int getK();

	public long getEntries();
	
	public void initSize(long size);
	
	public long getSize();

	public boolean put(CGATRingBuffer buffer, V value, boolean reverse);

	public boolean put(byte[] nseq, int start, V value, boolean reverse);

	public V get(CGATRingBuffer buffer, boolean reverse);

	public V get(byte[] nseq, int start, boolean reverse);

	public void visit(KMerStoreVisitor<V> visitor);
	
	public void optimize();
	
	public boolean isOptimized();
	
	public static <T extends Serializable> void save(KMerStore<T> store, File trieFile) throws IOException {
		ObjectOutputStream oOut = new ObjectOutputStream(StreamProvider.getOutputStreamForFile(trieFile));
		oOut.writeObject(store);
		oOut.close();		
	}
	
	public static <T extends Serializable> KMerStore<T> load(File filterFile)
			throws IOException, ClassNotFoundException {
		ObjectInputStream oOut = new ObjectInputStream(StreamProvider.getInputStreamForFile(filterFile));
		@SuppressWarnings("unchecked")
		KMerStore<T> res = (KMerStore<T>) oOut.readObject();
		oOut.close();
		return res;
	}
	
	public interface KMerStoreVisitor<V extends Serializable> {
		public void nextValue(KMerStore<V> trie, long kmer, V value);
	}
}
