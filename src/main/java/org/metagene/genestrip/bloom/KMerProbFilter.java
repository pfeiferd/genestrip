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

public interface KMerProbFilter extends Serializable {
    public void putLong(final long data);
    public boolean containsLong(final long data);
    public long ensureExpectedSize(long expectedInsertions, boolean enforceLarge);
    public void clear();
    public long getEntries();

    default public void save(File filterFile) throws IOException {
        try (ObjectOutputStream oOut = new ObjectOutputStream(StreamProvider.getOutputStreamForFile(filterFile))) {
            oOut.writeObject(this);
        }
    }

    public static KMerProbFilter load(InputStream is) throws IOException, ClassNotFoundException {
        try (ObjectInputStream oOut = new ObjectInputStream(is)) {
            return (KMerProbFilter) oOut.readObject();
        }
    }

    public static KMerProbFilter load(File filterFile) throws IOException, ClassNotFoundException {
        try (InputStream is = StreamProvider.getInputStreamForFile(filterFile)) {
            return load(is);
        }
    }
}