package org.metagene.genestrip.trie;

import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import org.metagene.genestrip.util.CGATRingBuffer;
import org.metagene.genestrip.util.StreamProvider;

public interface KMerStore<V extends Serializable> {
	public int getLen();

	public long getEntries();

	public void put(CGATRingBuffer buffer, V value, boolean reverse);

	public void put(byte[] nseq, int start, V value, boolean reverse);

	public V get(CGATRingBuffer buffer, boolean reverse);

	public V get(byte[] nseq, int start, boolean reverse);

	public void visit(KMerStoreVisitor<V> visitor, boolean reverse);
	
	public void optimize();
	
	public boolean isOptimized();

	public interface KMerStoreVisitor<V extends Serializable> {
		public void nextValue(KMerStore<V> trie, byte[] kmer, V value);
	}

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
}
