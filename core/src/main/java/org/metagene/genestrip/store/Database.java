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

/**
 * A Genestrip database: a k-mer {@link KMerStore} keyed by taxonomic id, the {@link SmallTaxTree} it
 * was built against, and the configuration properties describing how it was created. Persisted as a
 * ZIP holding the serialized database, its probabilistic pre-filter and a config-info entry (which
 * carries an MD5 fingerprint of the serialized database).
 */
public class Database implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * ZIP entry name of the serialized database.
     */
    public static final String DB_FILE = "db.ser";
    /**
     * ZIP entry name of the serialized probabilistic pre-filter.
     */
    public static final String INDEX_FILE = "bloom.ser";
    /**
     * ZIP entry name of the config-info properties.
     */
    public static final String CONFIG_INFO_FILE = "configInfo.properties";

    /**
     * The taxonomy tree this database was built against.
     */
    private final SmallTaxTree taxTree;
    /**
     * The k-mer store keyed by taxonomic id.
     */
    private final KMerStore<String> kmerStore;
    /**
     * Configuration properties describing how this database was created.
     */
    private Properties configInfo;

    /**
     * Creates a database from the given k-mer store, taxonomy tree and config-info properties.
     *
     * @param kmerStore  the k-mer store keyed by taxonomic id
     * @param taxTree    the taxonomy tree the store was built against
     * @param configInfo the configuration properties (a new empty {@link Properties} is used if {@code null})
     */
    public Database(KMerStore<String> kmerStore, SmallTaxTree taxTree, Properties configInfo) {
        this.kmerStore = kmerStore;
        this.taxTree = taxTree;
        this.configInfo = configInfo == null ? new Properties() : configInfo;
    }

    /**
     * Returns the k-mer store backing this database.
     *
     * @return the k-mer store keyed by taxonomic id
     */
    public KMerStore<String> getKmerStore() {
        return kmerStore;
    }

    /**
     * Assigns a store value index to every taxon in the taxonomy tree (see
     * {@link KMerStore#getAddValueIndex}).
     */
    public void initStoreIndices() {
        initStoreIndex(taxTree.getRoot());
    }

    /**
     * Recursively assigns store value indices to the given tax node and its descendants.
     *
     * @param node the tax node to start from (ignored if {@code null} or without a taxid)
     */
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

    /**
     * Returns a copy of the k-mer store whose taxid values are replaced by the corresponding
     * {@link SmallTaxIdNode}s of this database's tax tree.
     *
     * @return a converted copy of the k-mer store with tax node values
     */
    public KMerStore<SmallTaxIdNode> convertKMerStore() {
        return kmerStore.convertValues(new ValueConverter<String, SmallTaxIdNode>() {
            @Override
            public SmallTaxIdNode convertValue(String value) {
                return taxTree.getNodeByTaxId(value);
            }
        });
    }

    /**
     * Returns the taxonomy tree this database was built against.
     *
     * @return the taxonomy tree this database was built against
     */
    public SmallTaxTree getTaxTree() {
        return taxTree;
    }

    /**
     * Returns per-taxid store statistics.
     *
     * @return the number of stored k-mers per taxid, see {@link KMerStore#getNKmersPerTaxid()}.
     */
    public Object2LongMap<String> getStats() {
        return kmerStore.getNKmersPerTaxid();
    }

    /**
     * Saves this database (and its pre-filter, if any) to the given file, see {@link #save(OutputStream)}.
     *
     * @param file the destination file
     * @throws java.io.IOException if writing fails
     */
    public void save(File file) throws IOException {
        try (OutputStream out = new FileOutputStream(file)) {
            save(out);
        }
    }

    /**
     * Returns the configuration properties describing how this database was created.
     *
     * @return the configuration properties describing how this database was created
     */
    public Properties getConfigInfo() {
        return configInfo;
    }

    /**
     * Sets the configuration properties for this database.
     *
     * @param versionInfo the configuration properties to set
     */
    protected void setConfigInfo(Properties versionInfo) {
        this.configInfo = versionInfo;
    }

    /**
     * Writes this database to the stream as a ZIP with three entries: the serialized database, the
     * probabilistic pre-filter (only for a {@link TunableKMerStore}, else {@code null}) and the
     * config-info properties, the latter carrying an MD5 fingerprint of the serialized database bytes.
     *
     * @param os the destination stream
     * @throws java.io.IOException if writing fails
     */
    public void save(OutputStream os) throws IOException {
        try (ZipOutputStream zipOut = new ZipOutputStream(os)) {
            ZipEntry zipEntry;

            MessageDigest messageDigest = MessageDigest.getInstance("MD5");
            DigestOutputStream digo = new DigestOutputStream(zipOut, messageDigest);
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

            // We must NOT close digo, because it would close zipOut too.
            // We just leave digo hanging, it is just a FilterOutputStream - no worries if never closed.

            // From here on, we use zipOut directly, because only the core DB data shall be in the MD5 fingerprint.
            Properties configInfo = getConfigInfo();
            String md5 = Hex.encodeHexString(messageDigest.digest());
            configInfo.setProperty(GSProject.DB_MD5, md5);
            zipEntry = new ZipEntry(CONFIG_INFO_FILE);
            zipOut.putNextEntry(zipEntry);
            configInfo.store(zipOut, "Genestrip database configuration information");
            zipOut.closeEntry();
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Loads a database from the given file, see {@link #load(InputStream, boolean)}.
     *
     * @param file       the database file to read
     * @param withFilter whether to also load the probabilistic pre-filter
     * @return the loaded database
     * @throws java.io.IOException              if reading fails
     * @throws java.lang.ClassNotFoundException if a serialized class cannot be resolved
     */
    public static Database load(File file, boolean withFilter) throws IOException, ClassNotFoundException {
        try (InputStream is = new FileInputStream(file)) {
            return load(is, withFilter);
        }
    }

    /**
     * Loads a database from the given ZIP stream, restoring its store value indices. When
     * {@code withFilter} is set and the store is tunable, its probabilistic pre-filter is loaded too.
     *
     * @param is         the ZIP stream to read
     * @param withFilter whether to also load the probabilistic pre-filter
     * @return the loaded database
     * @throws java.lang.ClassNotFoundException if a serialized class cannot be resolved
     * @throws InvalidDatabaseClassException    if the serialized classes are incompatible with the current
     *                                          runtime (carrying the loaded config info for diagnostics).
     */
    public static Database load(InputStream is, boolean withFilter) throws IOException, ClassNotFoundException {
        Properties configInfo = new Properties();
        // The config entry (which holds the DB version/title) is stored AFTER the DB and index
        // entries (it carries the MD5 of the already-written DB). So a class-incompatibility error
        // while deserializing the DB or index is deferred until the loop has read the config -
        // otherwise InvalidDatabaseClassException would report an empty config ("unknown" version).
        InvalidClassException deferred = null;
        try {
            Database database = null;
            KMerProbFilter filter = null;
            try (ZipInputStream zis = new ZipInputStream(is)) {
                for (ZipEntry zipEntry = zis.getNextEntry(); zipEntry != null; zipEntry = zis.getNextEntry()) {
                    String entryName = zipEntry.getName();
                    if (entryName.equals(CONFIG_INFO_FILE)) {
                        configInfo.load(zis);
                        zis.closeEntry();
                    } else if (entryName.equals(DB_FILE) && database == null && deferred == null) {
                        ObjectInputStream oOut = new ObjectInputStream(zis);
                        try {
                            database = (Database) oOut.readObject();
                        } catch (InvalidClassException e) {
                            deferred = e;
                        }
                        zis.closeEntry();
                    } else if (entryName.equals(INDEX_FILE) && withFilter && filter == null && deferred == null) {
                        ObjectInputStream oOut = new ObjectInputStream(zis);
                        try {
                            filter = (KMerProbFilter) oOut.readObject();
                        } catch (InvalidClassException e) {
                            deferred = e;
                        }
                        zis.closeEntry();
                    }
                }
                if (deferred != null) {
                    // configInfo is now populated; re-throw the raw exception and let the outer
                    // catch wrap it (once) into an InvalidDatabaseClassException carrying the config.
                    throw deferred;
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

    /**
     * Loads only the config-info properties from the given database file.
     *
     * @param file the database file to read
     * @return the config-info properties
     * @throws java.io.IOException if reading fails
     */
    public static Properties loadConfigInfo(File file) throws IOException {
        try (InputStream is = new FileInputStream(file)) {
            return loadConfigInfo(is);
        }
    }

    /**
     * Reads only the config-info entry from the given database ZIP stream, without deserializing the
     * database itself.
     *
     * @param is the ZIP stream to read
     * @return the config-info properties
     * @throws java.io.IOException if reading fails
     */
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
