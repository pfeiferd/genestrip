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
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.tax.TaxIdCollector;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.StreamProvider;

public class KMerStoreWrapper implements Serializable {
	private static final long serialVersionUID = 1L;

	private final KMerSortedArray<String> kmerStore;
	private final List<TaxIdNode> taxids;
	private final Map<String, StoreStatsPerTaxid> storeStats;

	public KMerStoreWrapper(KMerSortedArray<String> kmerStore, Set<TaxIdNode> taxids,
			Map<String, StoreStatsPerTaxid> storeStats) {
		this(kmerStore, TaxIdCollector.nodesAsShallowCopies(TaxIdCollector.sortNodes(taxids)), storeStats);
	}

	public KMerStoreWrapper(KMerSortedArray<String> kmerStore, List<TaxIdNode> taxids,
			Map<String, StoreStatsPerTaxid> storeStats) {
		this.kmerStore = kmerStore;
		this.taxids = taxids;
		this.storeStats = storeStats;
	}

	public KMerSortedArray<String> getKmerStore() {
		return kmerStore;
	}

	public List<TaxIdNode> getTaxids() {
		return taxids;
	}

	public Map<String, StoreStatsPerTaxid> getStoreStats() {
		return storeStats;
	}

	public void save(File file) throws IOException {
		ObjectOutputStream oOut = new ObjectOutputStream(StreamProvider.getOutputStreamForFile(file));
		oOut.writeObject(this);
		oOut.close();
	}

	public static KMerStoreWrapper load(File file) throws IOException, ClassNotFoundException {
		ObjectInputStream oOut = new ObjectInputStream(StreamProvider.getInputStreamForFile(file));
		KMerStoreWrapper res = (KMerStoreWrapper) oOut.readObject();
		oOut.close();
		return res;
	}

	public static class StoreStatsPerTaxid implements Serializable {
		private static final long serialVersionUID = 1L;
		
		public long totalKMers;
		public long storedKMers;
		public long contigs;
		public long maxContigLen;

		public long getContigs() {
			return contigs;
		}

		public long getMaxContigLen() {
			return maxContigLen;
		}

		public long getStoredKMers() {
			return storedKMers;
		}

		public long getTotalKMers() {
			return totalKMers;
		}
	}
}
