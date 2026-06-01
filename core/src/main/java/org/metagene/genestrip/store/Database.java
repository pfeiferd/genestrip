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
import java.security.DigestInputStream;
import java.security.DigestOutputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Properties;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

import org.metagene.genestrip.GSProject;
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
    public static final String MD5_FILE = "md5";

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

    public void initStoreIndices() {
        initStoreIndex(taxTree.getRoot());
    }

    protected void initStoreIndex(SmallTaxIdNode node) {
        if (node == null || node.getTaxId() == null) {
            return;
        }
        int index = kmerStore.getAddValueIndex(node.getTaxId());
        node.setStoreIndex(index);
        SmallTaxIdNode[] subnodes = node.getSubNodes();
        if (subnodes != null) {
            for (int i = 0; i < subnodes.length; i++) {
                initStoreIndex(subnodes[i]);
            }
        }
    }

    public void ensureAllTreeNodesInDB() {
        ensureAllTreeNodesInDB(taxTree.getRoot());
    }

    protected void ensureAllTreeNodesInDB(SmallTaxIdNode node) {
        if (node == null || node.getTaxId() == null) {
            return;
        }
        int index = kmerStore.getAddValueIndex(node.getTaxId());
        node.setStoreIndex(index);
        SmallTaxIdNode[] subnodes = node.getSubNodes();
        if (subnodes != null) {
            for (int i = 0; i < subnodes.length; i++) {
                ensureAllTreeNodesInDB(subnodes[i]);
            }
        }
    }

    public KMerSortedArray<SmallTaxIdNode> convertKMerStore() {
        return new KMerSortedArray<>(kmerStore, new ValueConverter<String, SmallTaxIdNode>() {
            @Override
            public SmallTaxIdNode convertValue(String value) {
                return taxTree.getNodeByTaxId(value);
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
        DigestOutputStream digo = null;
        try (ZipOutputStream zipOut = new ZipOutputStream(os)) {
            ZipEntry zipEntry = null;
            Properties configInfo = getConfigInfo();
            if (configInfo != null) {
                zipEntry = new ZipEntry(CONFIG_INFO_FILE);
                zipOut.putNextEntry(zipEntry);
                // Zip out here, cause only DB itself shall be in the md5 finger print.
                configInfo.store(zipOut, "Genestrip database configuration information");
                zipOut.closeEntry();
            }

            MessageDigest messageDigest = MessageDigest.getInstance("MD5");
            // No try() here cause that would close the zip stream too.
            digo = new DigestOutputStream(zipOut, messageDigest);
            zipEntry = new ZipEntry(DB_FILE);
            zipOut.putNextEntry(zipEntry);
            ObjectOutputStream oOut = new ObjectOutputStream(digo);
            oOut.writeObject(this);
            zipOut.closeEntry();

            zipEntry = new ZipEntry(INDEX_FILE);
            zipOut.putNextEntry(zipEntry);
            oOut = new ObjectOutputStream(digo);
            oOut.writeObject(kmerStore.getFilter());
            digo.flush(); // Make sure it all got written.
            zipOut.closeEntry();

            String md5 = messageDigest.digest().toString();
            zipEntry = new ZipEntry(MD5_FILE);
            zipOut.putNextEntry(zipEntry);
            zipOut.write(md5.getBytes());
            zipOut.closeEntry();
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException(e);
        } finally {
            if (digo != null) {
                digo.close();
            }
        }
    }

    public static Database load(File file, boolean withFilter) throws IOException, ClassNotFoundException {
        try (InputStream is = new FileInputStream(file)) {
            return load(is, withFilter);
        }
    }

    public static Database load(InputStream is, boolean withFilter) throws IOException, ClassNotFoundException {
        Properties configInfo = new Properties();
        Database database = null;
        KMerProbFilter filter = null;
        try (ZipInputStream zis = new ZipInputStream(is)) {
            for (ZipEntry zipEntry = zis.getNextEntry(); zipEntry != null; zipEntry = zis.getNextEntry()) {
                String entryName = zipEntry.getName();
                // We want the version info as the first thing.
                // We assume it got stored as first thing and zip maintains this order in the
                // underlying zip file.
                if (entryName.equals(CONFIG_INFO_FILE)) {
                    configInfo.load(zis);
                    zis.closeEntry();
                } else if (entryName.equals(DB_FILE) && database == null) {
                    ObjectInputStream oOut = new ObjectInputStream(zis);
                    database = (Database) oOut.readObject();
                    zis.closeEntry();
                } else if (entryName.equals(INDEX_FILE) && withFilter && filter == null) {
                    ObjectInputStream oOut = new ObjectInputStream(zis);
                    filter = (KMerProbFilter) oOut.readObject();
                    zis.closeEntry();
                } else if (entryName.equals(MD5_FILE)) {
                    byte[] md5b = new byte[64]; // Large enough buffer for sure
                    int size = zis.read(md5b);
                    zis.closeEntry();
                    String md5 = new String(md5b, 0, size);
                    configInfo.setProperty(GSProject.DB_MD5, md5);
                }
            }
            database.setConfigInfo(configInfo);
            database.initStoreIndices();
            if (filter != null) {
                database.getKmerStore().setFilter(filter);
            }
        }
        return database;
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
