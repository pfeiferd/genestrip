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
     * Adds the given k-mer to the filter unless it is (probably) already present, and reports whether
     * it was newly added. This is required to be equivalent to a {@code if (!containsLong(data))
     * putLong(data)} sequence, but implementations combine it into a single pass and, where they
     * support it, make it safe for concurrent (lock-free) use. Deliberately left without a default:
     * a naive {@code containsLong}/{@code putLong} fallback would not be thread-safe, so every
     * implementation must decide and document its own concurrency behaviour rather than silently
     * inherit an unsafe one (see {@link AbstractKMerBloomFilter} and {@link BlockedKMerBloomFilter}).
     * <p>
     * For the concurrent implementations, two threads inserting the same absent k-mer may both
     * observe it as new, so under concurrent use the returned flag may marginally over-count the
     * newly-added k-mers; membership answers are never affected.
     *
     * @param data the k-mer, encoded as a {@code long}, to add
     * @return {@code true} if the k-mer was not already present, {@code false} otherwise
     */
    public boolean putLong(final long data);

    /**
     * Tests whether the given k-mer is (probably) contained in the filter.
     *
     * @param data the k-mer, encoded as a {@code long}, to test
     * @return {@code false} if the k-mer is definitely absent, {@code true} if it is (probably) present.
     */
    public boolean containsLong(final long data);
    /**
     * Removes all entries, keeping the current capacity.
     */
    public void clear();

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