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

import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.concurrent.atomic.LongAdder;

/**
 * Base class for {@link KMerProbFilter} implementations that support a concurrent {@link
 * #putLongIfAbsent(long)} insert path. It factors out the (subtle) lifecycle of the concurrent
 * entry counter that path relies on: a {@link LongAdder} accumulating the newly-added k-mers without
 * contending on a single shared counter field.
 * <p>
 * The plain {@code entries} counter (written by {@code putLong}) is deliberately <em>not</em> pulled
 * up here: each subclass keeps declaring it itself so that the on-disk serialization form of existing
 * filters is preserved. Subclasses combine their {@code entries} with {@link #concurrentEntryCount()}
 * in {@code getEntries()} and fold {@link #drainConcurrentEntries()} into {@code entries} in their
 * {@code writeObject}. The adder itself is transient, so adding this class to the hierarchy does not
 * change the serialized form.
 */
public abstract class AbstractCountingKMerProbFilter implements KMerProbFilter {
	private static final long serialVersionUID = 1L;

	/**
	 * Concurrent counter of the new k-mers added via {@link #putLongIfAbsent(long)}. Kept separate
	 * from the subclass's plain {@code entries} counter so the lock-free insert path does not need to
	 * synchronize on a shared field; subclasses fold it into {@code entries} on serialization.
	 * Transient because it only matters while the filter is being populated.
	 */
	protected transient LongAdder concurrentEntries;

	/**
	 * Creates the filter with a fresh, zeroed concurrent counter.
	 */
	protected AbstractCountingKMerProbFilter() {
		concurrentEntries = new LongAdder();
	}

	/**
	 * Restores the transient concurrent counter after deserialization, so a reloaded filter can be
	 * populated again via {@link #putLongIfAbsent(long)} without a {@code NullPointerException}.
	 *
	 * @param in the stream to read from
	 * @throws IOException            if reading fails
	 * @throws ClassNotFoundException if a serialized class cannot be resolved
	 */
	private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException {
		in.defaultReadObject();
		ensureConcurrentEntries();
	}

	/**
	 * Ensures the concurrent counter exists, recreating it if it was dropped (e.g. after
	 * deserialization). Subclasses call this before populating the filter.
	 */
	protected final void ensureConcurrentEntries() {
		if (concurrentEntries == null) {
			concurrentEntries = new LongAdder();
		}
	}

	/**
	 * Resets the concurrent counter to zero (recreating it if necessary). Subclasses call this when
	 * (re)sizing a filter that should start counting from scratch.
	 */
	protected final void resetConcurrentEntries() {
		if (concurrentEntries == null) {
			concurrentEntries = new LongAdder();
		} else {
			concurrentEntries.reset();
		}
	}

	/**
	 * Records that one k-mer was newly added via the lock-free insert path. Safe to call concurrently.
	 */
	protected final void countAdded() {
		concurrentEntries.increment();
	}

	/**
	 * Returns the number of k-mers counted via the concurrent insert path so far, for a subclass to
	 * add to its plain {@code entries} counter in {@code getEntries()}.
	 *
	 * @return the current concurrent entry count, or {@code 0} if no counter is present
	 */
	protected final long concurrentEntryCount() {
		return concurrentEntries == null ? 0 : concurrentEntries.sum();
	}

	/**
	 * Returns the number of k-mers counted via the concurrent insert path and resets that counter to
	 * zero, for a subclass to fold into its plain {@code entries} counter in {@code writeObject}
	 * before serialization.
	 *
	 * @return the concurrent entry count that was drained, or {@code 0} if no counter is present
	 */
	protected final long drainConcurrentEntries() {
		return concurrentEntries == null ? 0 : concurrentEntries.sumThenReset();
	}
}
