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

import org.apache.commons.codec.binary.Hex;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.bloom.KMerProbFilter;
import org.metagene.genestrip.store.KMerStore.ValueConverter;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class Database implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String DB_FILE = "db.ser";
    public static final String INDEX_FILE = "bloom.ser";
    public static final String CONFIG_INFO_FILE = "configInfo.properties";

    private final SmallTaxTree taxTree;
    private final KMerStore<String> kmerStore;
    private Properties configInfo;

    public Database(KMerStore<String> kmerStore, SmallTaxTree taxTree, Properties configInfo) {
        this.kmerStore = kmerStore;
        this.taxTree = taxTree;
        this.configInfo = configInfo == null ? new Properties() : configInfo;
    }

    public KMerStore<String> getKmerStore() {
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

    public KMerStore<SmallTaxIdNode> convertKMerStore() {
        return kmerStore.convertValues(new ValueConverter<String, SmallTaxIdNode>() {
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

            MessageDigest messageDigest = MessageDigest.getInstance("MD5");
            // No try() here cause that would close the zip stream too.
            {
                digo = new DigestOutputStream(zipOut, messageDigest);
                zipEntry = new ZipEntry(DB_FILE);
                zipOut.putNextEntry(zipEntry);
                ObjectOutputStream oOut = new ObjectOutputStream(digo);
                oOut.writeObject(this);
                zipOut.closeEntry();

                zipEntry = new ZipEntry(INDEX_FILE);
                zipOut.putNextEntry(zipEntry);
                oOut = new ObjectOutputStream(digo);
                // Only tunable stores expose a probabilistic pre-filter; for others none is written.
                KMerProbFilter filter = kmerStore instanceof TunableKMerStore
                        ? ((TunableKMerStore<String>) kmerStore).getFilter() : null;
                oOut.writeObject(filter);
                digo.flush(); // Make sure it all got written.
                zipOut.closeEntry();
            }

            // Zip out here, cause only DB itself shall be in the md5 finger print.
            Properties configInfo = getConfigInfo();
            String md5 = Hex.encodeHexString(messageDigest.digest());
            configInfo.setProperty(GSProject.DB_MD5, md5);
            zipEntry = new ZipEntry(CONFIG_INFO_FILE);
            zipOut.putNextEntry(zipEntry);
            configInfo.store(zipOut, "Genestrip database configuration information");
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
        try {
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
                    }
                }
                database.setConfigInfo(configInfo);
                database.initStoreIndices();
                if (filter != null && database.getKmerStore() instanceof TunableKMerStore) {
                    ((TunableKMerStore<String>) database.getKmerStore()).setFilter(filter);
                }
            }
            return database;
        } catch (InvalidClassException e) {
            throw new InvalidDatabaseClassException(e.classname, e.getMessage(), configInfo, e);
        }
    }

    public static Properties loadConfigInfo(File file) throws IOException {
        try (InputStream is = new FileInputStream(file)) {
            return loadConfigInfo(is);
        }
    }

    public static Properties loadConfigInfo(InputStream is) throws IOException {
        Properties configInfo = new Properties();
        try (ZipInputStream zis = new ZipInputStream(is)) {
            for (ZipEntry zipEntry = zis.getNextEntry(); zipEntry != null; zipEntry = zis.getNextEntry()) {
                if (zipEntry.getName().equals(CONFIG_INFO_FILE)) {
                    configInfo.load(zis);
                    zis.closeEntry();
                }
            }
        }
        return configInfo;
    }
}
