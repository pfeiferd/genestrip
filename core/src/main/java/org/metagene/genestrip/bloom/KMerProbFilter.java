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
package org.metagene.genestrip.bloom;

import org.metagene.genestrip.io.StreamProvider;

import java.io.*;

/**
 * A probabilistic set of k-mers (each encoded as a {@code long}), i.e. a Bloom-style filter used to
 * pre-screen k-mer lookups: {@link #containsLong} may return a false positive but never a false
 * negative.
 */
public interface KMerProbFilter extends Serializable {
    /**
     * Adds the given k-mer to the filter.
     *
     * @param data the k-mer, encoded as a {@code long}, to add
     */
    public void putLong(final long data);
    /**
     * Tests whether the given k-mer is (probably) contained in the filter.
     *
     * @param data the k-mer, encoded as a {@code long}, to test
     * @return {@code false} if the k-mer is definitely absent, {@code true} if it is (probably) present.
     */
    public boolean containsLong(final long data);
    /**
     * (Re)sizes the filter for the given expected number of insertions, using large backing storage
     * when required or when {@code enforceLarge} is set.
     *
     * @param expectedInsertions the expected number of k-mers to be inserted
     * @param enforceLarge        whether to force the use of large backing storage
     * @return the resulting number of bits.
     */
    public long ensureExpectedSize(long expectedInsertions, boolean enforceLarge);
    /**
     * Removes all entries, keeping the current capacity.
     */
    public void clear();
    /**
     * Returns the number of k-mers that have been added to the filter.
     *
     * @return the number of inserted entries
     */
    public long getEntries();

    /**
     * Serializes this filter to the given file.
     *
     * @param filterFile the file to write the serialized filter to
     * @throws IOException if the filter cannot be written
     */
    default public void save(File filterFile) throws IOException {
        try (ObjectOutputStream oOut = new ObjectOutputStream(StreamProvider.getOutputStreamForFile(filterFile))) {
            oOut.writeObject(this);
        }
    }

    /**
     * Deserializes a filter from the given stream.
     *
     * @param is the stream to read the serialized filter from
     * @return the deserialized filter
     * @throws IOException            if the stream cannot be read
     * @throws ClassNotFoundException if the serialized filter class cannot be resolved
     */
    public static KMerProbFilter load(InputStream is) throws IOException, ClassNotFoundException {
        try (ObjectInputStream oOut = new ObjectInputStream(is)) {
            return (KMerProbFilter) oOut.readObject();
        }
    }

    /**
     * Deserializes a filter from the given file.
     *
     * @param filterFile the file to read the serialized filter from
     * @return the deserialized filter
     * @throws IOException            if the file cannot be read
     * @throws ClassNotFoundException if the serialized filter class cannot be resolved
     */
    public static KMerProbFilter load(File filterFile) throws IOException, ClassNotFoundException {
        try (InputStream is = StreamProvider.getInputStreamForFile(filterFile)) {
            return load(is);
        }
    }
}