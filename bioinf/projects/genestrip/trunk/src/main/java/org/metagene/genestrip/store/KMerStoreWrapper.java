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

import org.metagene.genestrip.bloom.MurmurCGATBloomFilter;
import org.metagene.genestrip.io.StreamProvider;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class KMerStoreWrapper implements Serializable {
	private static final long serialVersionUID = 1L;

	private final KMerSortedArray<String> kmerStore;

	public KMerStoreWrapper(KMerSortedArray<String> kmerStore) {
		this.kmerStore = kmerStore;
	}

	public KMerSortedArray<String> getKmerStore() {
		return kmerStore;
	}
	
	public Object2LongMap<String> getStats() {
		return kmerStore.getStats();
	}

	public void save(File file, File filterFile) throws IOException {
		ObjectOutputStream oOut = new ObjectOutputStream(StreamProvider.getOutputStreamForFile(file));
		oOut.writeObject(this);
		oOut.close();
		if (filterFile != null) {
			kmerStore.getFilter().save(filterFile);
		}
	}

	public static KMerStoreWrapper load(File file, File filterFile) throws IOException, ClassNotFoundException {
		ObjectInputStream oOut = new ObjectInputStream(StreamProvider.getInputStreamForFile(file));
		KMerStoreWrapper res = (KMerStoreWrapper) oOut.readObject();
		oOut.close();
		if (filterFile != null) {
			MurmurCGATBloomFilter filter = MurmurCGATBloomFilter.load(filterFile);
			res.getKmerStore().setFilter(filter);
		}
		
		return res;
	}
}
