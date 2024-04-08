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
package org.metagene.genestrip.goals;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import org.metagene.genestrip.GSConfig;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.FastqKMerMatcher2;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.match.ResultReporter;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.store.KMerStoreWrapper;
import org.metagene.genestrip.store.KMerUniqueCounterBits;
import org.metagene.genestrip.store.KMerSortedArray.ValueConverter;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ExecutorServiceBundle;

public class MatchGoal extends FileListGoal<GSProject> {
	private final File fastq;
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final ObjectGoal<KMerStoreWrapper, GSProject> storeGoal;
	private final boolean writedFiltered;
	private final ExecutorServiceBundle bundle;
	private FastqKMerMatcher2 matcher;

	@SafeVarargs
	public MatchGoal(GSProject project, String name, File fastq, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			ObjectGoal<KMerStoreWrapper, GSProject> storeGoal, boolean writeFiltered, ExecutorServiceBundle bundle,
			Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, null, fastq.getName(), FileType.CSV, false),
				Goal.append(deps, taxTreeGoal, storeGoal));
		this.fastq = fastq;
		this.taxTreeGoal = taxTreeGoal;
		this.storeGoal = storeGoal;
		this.writedFiltered = writeFiltered;
		this.bundle = bundle;
	}

	@Override
	protected void makeFile(File file) {
		matcher = null;
		try {
			File filteredFile = null;
			File krakenOutStyleFile = null;
			if (writedFiltered) {
				filteredFile = getProject().getOutputFile(getName(), null, fastq.getName(), FileType.FASTQ_RES, true);
				krakenOutStyleFile = getProject().getOutputFile(getName(), null, fastq.getName(),
						FileType.KRAKEN_OUT_RES, false);
			}

			KMerStoreWrapper wrapper = storeGoal.get();
			wrapper.getKmerStore().setUseFilter(getProject().isUseBloomFilterForMatch());

			GSConfig config = getProject().getConfig();

			KMerSortedArray<TaxIdNode> store = new KMerSortedArray<TaxIdNode>(wrapper.getKmerStore(),
					new ValueConverter<String, TaxIdNode>() {
						@Override
						public TaxIdNode convertValue(String value) {
							return taxTreeGoal.get().getNodeByTaxId(value);
						}
					});

			matcher = createMatcher(store, taxTreeGoal.get(), bundle);
			MatchingResult res = matcher.runMatcher(new StreamingFileResource(fastq), filteredFile, krakenOutStyleFile,
					config.isCountUniqueKMers()
							? new KMerUniqueCounterBits(wrapper.getKmerStore(), config.isMatchWithKMerCounts())
							: null);
			matcher.dump();

			writeOutputFile(file, res, wrapper);
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			dumpMatcher();
		}
	}

	@Override
	public void dump() {
		super.dump();
		dumpMatcher();
	}

	protected void dumpMatcher() {
		if (matcher != null) {
			matcher.dump();
			matcher = null;
		}
	}

	protected void writeOutputFile(File file, MatchingResult result, KMerStoreWrapper wrapper) throws IOException {
		PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file));
		new ResultReporter(taxTreeGoal.get(), getProject().getConfig().getNormalizedKMersFactor())
				.printMatchResult(result, out, wrapper);
		out.close();
	}

	protected FastqKMerMatcher2 createMatcher(KMerSortedArray<TaxIdNode> store, TaxTree taxTree,
			ExecutorServiceBundle bundle) {
		GSConfig config = getProject().getConfig();

		return new FastqKMerMatcher2(store, config.getMaxReadSizeBytes(), config.getThreadQueueSize(), bundle,
				config.getMaxKMerResCounts(), getProject().isClassifyReads() ? taxTree : null,
				config.getMaxClassificationPaths(), getProject().getMaxReadTaxErrorCount());
	}
}
