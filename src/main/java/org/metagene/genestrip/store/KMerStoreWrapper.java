package org.metagene.genestrip.store;

import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.List;
import java.util.Set;

import org.metagene.genestrip.tax.TaxIdCollector;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.StreamProvider;

public class KMerStoreWrapper implements Serializable {
	private static final long serialVersionUID = 1L;
	
	private final KMerStore<String> kmerStore;
	private final List<TaxIdNode> taxids;
	
	public KMerStoreWrapper(KMerStore<String> kmerStore, Set<TaxIdNode> taxids) {
		this(kmerStore, TaxIdCollector.nodesAsShallowCopies(TaxIdCollector.sortNodes(taxids)));
	}
	
	public KMerStoreWrapper(KMerStore<String> kmerStore, List<TaxIdNode> taxids) {
		this.kmerStore = kmerStore;
		this.taxids = taxids;		
	}
	
	public KMerStore<String> getKmerStore() {
		return kmerStore;
	}
	
	public List<TaxIdNode> getTaxids() {
		return taxids;
	}
	
	public void save(File file) throws IOException {
		ObjectOutputStream oOut = new ObjectOutputStream(StreamProvider.getOutputStreamForFile(file));
		oOut.writeObject(this);
		oOut.close();		
	}
	
	public static KMerStoreWrapper load(File file)
			throws IOException, ClassNotFoundException {
		ObjectInputStream oOut = new ObjectInputStream(StreamProvider.getInputStreamForFile(file));
		KMerStoreWrapper res = (KMerStoreWrapper) oOut.readObject();
		oOut.close();
		return res;
	}
}
