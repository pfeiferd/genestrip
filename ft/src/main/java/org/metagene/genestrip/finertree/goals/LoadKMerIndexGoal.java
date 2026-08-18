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
package org.metagene.genestrip.finertree.goals;

import me.tongfei.progressbar.ProgressBar;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.probfilter.ProbFilter;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InvalidClassException;

/**
 * Provides the {@link ProbFilter} k-mer index, either by reusing the in-memory filter
 * already produced by the bloom goal or, if that goal was not made, by deserializing the filter from
 * the file written by {@link StoreKMerIndexGoal}.
 *
 * @param <P> the concrete {@link FTProject} type this goal operates on
 */
public class LoadKMerIndexGoal<P extends FTProject> extends ObjectGoal<ProbFilter, P> implements Goal.LogHeapInfo {
	private final ObjectGoal<ProbFilter, P> bloomIndex;
	private final File dbFile;

	/**
	 * Creates the goal that loads or reuses the k-mer index bloom filter.
	 *
	 * @param project        the project this goal belongs to
	 * @param bloomIndex     goal supplying the freshly built in-memory filter, if available
	 * @param storeIndexGoal file goal supplying the serialized filter file to load otherwise
	 * @param dependencies   further goals this goal depends on
	 */
	@SafeVarargs
	public LoadKMerIndexGoal(P project, ObjectGoal<ProbFilter, P> bloomIndex,
                             FileGoal<P> storeIndexGoal, Goal<P>... dependencies) {
		super(project, FTGoalKey.LOAD_KMER_INDEX, append(dependencies, bloomIndex, storeIndexGoal));
		this.bloomIndex = bloomIndex;
		this.dbFile = storeIndexGoal.getFile();
	}

	/**
	 * Loads the k-mer index and publishes it as this goal's result. When the in-memory filter is not
	 * available and the progress bar is enabled, the file is read through a progress-tracking stream;
	 * otherwise it is read directly.
	 */
	@Override
	protected void doMakeThis() {
		try {
			if (booleanConfigValue(GSConfigKey.PROGRESS_BAR) && !bloomIndex.isMade()) {
				try (StreamingResource.StreamAccess sa = new StreamingFileResource(dbFile, false).openStream()) {
					try (ProgressBar pb = GSProgressBarCreator.newGSProgressBar(getKey().getName(), sa, null)) {
						doLoadIndex(sa.getInputStream());
					}
				}
			}
			else {
				try (InputStream is = StreamProvider.getInputStreamForFile(dbFile)) {
					doLoadIndex(is);
				}
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Sets this goal's result to the in-memory filter from the bloom goal if it was made, or else
	 * deserializes the filter from the given stream.
	 *
	 * @param stream the stream to deserialize the filter from when no in-memory filter exists
	 */
	protected void doLoadIndex(InputStream stream) {
		try {
			ProbFilter filter = bloomIndex.isMade() ? bloomIndex.get()
					: (ProbFilter) ProbFilter.load(stream);
			set(filter);
		} catch (ClassNotFoundException e) {
			throw new RuntimeException(e);
		} catch (InvalidClassException e) {
			throw new RuntimeException("Index file version does not match genestrip library version.", e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
