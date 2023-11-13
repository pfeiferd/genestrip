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
package org.metagene.genestrip.store;

import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.util.CGATRingBuffer;

public interface KMerStore<V extends Serializable> extends Serializable {
	public void initSize(long size);
	
	public int getK();

	public long getEntries();

	public long getSize();

	/**
	 * @deprecated
	 */
	public boolean put(CGATRingBuffer buffer, V value, boolean reverse);

	/**
	 * @deprecated
	 */
	public boolean put(byte[] nseq, int start, V value, boolean reverse);

	/**
	 * @deprecated
	 */
	public V get(CGATRingBuffer buffer, boolean reverse);

	/**
	 * @deprecated
	 */
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
