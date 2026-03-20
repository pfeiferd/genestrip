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

import java.io.*;
import java.text.DateFormat;
import java.util.Date;
import java.util.Locale;
import java.util.Properties;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

import org.metagene.genestrip.bloom.KMerProbFilter;
import org.metagene.genestrip.store.KMerSortedArray.ValueConverter;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class Database implements Serializable {
	private static final long serialVersionUID = 1L;

	public static final String DB_FILE = "db.ser";
	public static final String INDEX_FILE = "bloom.ser";
	public static final String CONFIG_INFO_FILE = "configInfo.properties";

	private final SmallTaxTree taxTree;
	private final KMerSortedArray<String> kmerStore;
	private Properties configInfo;

	public Database(KMerSortedArray<String> kmerStore, SmallTaxTree taxTree, Properties configInfo) {
		this.kmerStore = kmerStore;
		this.taxTree = taxTree;
		this.configInfo = configInfo;
	}

	public KMerSortedArray<String> getKmerStore() {
		return kmerStore;
	}

	public KMerSortedArray<SmallTaxIdNode> convertKMerStore() {
		return new KMerSortedArray<SmallTaxIdNode>(kmerStore, new ValueConverter<String, SmallTaxIdNode>() {
			@Override
			public SmallTaxIdNode convertValue(String value) {
				SmallTaxIdNode node = taxTree.getNodeByTaxId(value);
				if (node != null && node.getStoreIndex() == -1) {
					node.setStoreIndex(kmerStore.getIndexForValue(value));
				}
				return node;
			}
		});
	}

	public SmallTaxTree getTaxTree() {
		return taxTree;
	}

	public Object2LongMap<String> getStats() {
		return kmerStore.getFixedNKmersPerTaxid();
	}

	public void save(File file) throws IOException {
		try (OutputStream out = new FileOutputStream(file)) {
			save(out);
		}
	}

	public Properties getConfigInfo() {
		return configInfo;
	}

	protected void setConfigInfo(Properties versionInfo) {
		this.configInfo = versionInfo;
	}

	public void save(OutputStream os) throws IOException {
		try (ZipOutputStream zipOut = new ZipOutputStream(os)) {
			ZipEntry zipEntry = null;
			Properties versionInfo = getConfigInfo();
			if (versionInfo != null) {
				zipEntry = new ZipEntry(CONFIG_INFO_FILE);
				zipOut.putNextEntry(zipEntry);
				getConfigInfo().store(zipOut, "Genestrip database version information");
				zipOut.closeEntry();
			}

			zipEntry = new ZipEntry(DB_FILE);
			zipOut.putNextEntry(zipEntry);
			ObjectOutputStream oOut = new ObjectOutputStream(zipOut);
			oOut.writeObject(this);
			zipOut.closeEntry();

			zipEntry = new ZipEntry(INDEX_FILE);
			zipOut.putNextEntry(zipEntry);
			oOut = new ObjectOutputStream(zipOut);
			oOut.writeObject(kmerStore.getFilter());
			zipOut.closeEntry();
		}
	}

	public static Database load(File file, boolean withFilter) throws IOException, ClassNotFoundException {
		try (InputStream is = new FileInputStream(file)) {
			return load(is, withFilter);
		}
	}

	public static Database load(InputStream is, boolean withFilter) throws IOException, ClassNotFoundException {
		Properties configInfo = new Properties();
		Database res = null;
		KMerProbFilter filter = null;
		try {
			try (ZipInputStream zis = new ZipInputStream(is)) {
				for (ZipEntry zipEntry = zis.getNextEntry(); zipEntry != null; zipEntry = zis.getNextEntry()) {
					String entryName = zipEntry.getName();
					// We want the version info as the first thing.
					// We assume it got stored as first thing and zip maintains this order in the
					// underlying zip file.
					if (entryName.equals(CONFIG_INFO_FILE) && configInfo == null) {
						configInfo = new Properties();
						configInfo.load(zis);
						zis.closeEntry();
					} else if (entryName.equals(DB_FILE) && res == null) {
						ObjectInputStream oOut = new ObjectInputStream(zis);
						res = (Database) oOut.readObject();
						zis.closeEntry();
					} else if (entryName.equals(INDEX_FILE) && withFilter && filter == null) {
						ObjectInputStream oOut = new ObjectInputStream(zis);
						filter = (KMerProbFilter) oOut.readObject();
						zis.closeEntry();
					}
				}
				res.setConfigInfo(configInfo);
				if (filter != null) {
					res.getKmerStore().setFilter(filter);
				}
			}
		} catch (InvalidClassException e) {
			throw new InvalidDatabaseClassException(e.classname, e.getMessage(), configInfo);
		}
		return res;
	}

	/* Should never be used - so remove it (?)
	public static MurmurCGATBloomFilter loadFilter(File file) throws IOException, ClassNotFoundException {
		try (InputStream is = new FileInputStream(file)) {
			return loadFilter(is);
		}
	}

	public static MurmurCGATBloomFilter loadFilter(InputStream is) throws IOException, ClassNotFoundException {
		MurmurCGATBloomFilter filter = null;
		try (ZipInputStream zis = new ZipInputStream(is)) {
			for (ZipEntry zipEntry = zis.getNextEntry(); zipEntry != null; zipEntry = zis.getNextEntry()) {
				String entryName = zipEntry.getName();
				if (entryName.equals(INDEX_FILE) && filter == null) {
					ObjectInputStream oOut = new ObjectInputStream(zis);
					filter = (MurmurCGATBloomFilter) oOut.readObject();
					zis.closeEntry();
					break;
				}
			}
		}
		return filter;
	}
	 */
}
