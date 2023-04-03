package org.metagene.genestrip.goals;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.kraken.KrakenExecutor;
import org.metagene.genestrip.kraken.KrakenResultFastqMergeListener;
import org.metagene.genestrip.kraken.KrakenResultFastqMerger;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayToString;
import org.metagene.genestrip.util.CountingDigitTrie;

public class KrakenResCountGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal;

	@SafeVarargs
	public KrakenResCountGoal(GSProject project, ObjectGoal<Set<TaxIdNode>, GSProject> taxNodesGoal,
			Goal<GSProject>... dependencies) {
		super(project, "krakenrescount", project.getTaxCountsFile("krakenrescount"), dependencies);
		this.taxNodesGoal = taxNodesGoal;
	}

	@Override
	protected void makeFile(File file) throws IOException {
		CountingDigitTrie kmerCountTrie = new CountingDigitTrie();
		CountingDigitTrie classCountTrie = new CountingDigitTrie();

		final Set<String> taxIds;

		if (taxNodesGoal != null) {
			taxIds = new HashSet<String>();
			for (TaxIdNode node : taxNodesGoal.get()) {
				taxIds.add(node.getTaxId());
			}
		} else {
			taxIds = null;
		}

		KrakenExecutor krakenExecutor = new KrakenExecutor(getProject().getConfig().getKrakenBinFolder(),
				getProject().getConfig().getKrakenExecExpr()) {
			@Override
			protected void handleOutputStream(InputStream stream, File outFile) throws IOException {
				KrakenResultFastqMerger parser = new KrakenResultFastqMerger(
						getProject().getConfig().getMaxReadSizeBytes());

				parser.process(new BufferedInputStream(stream), null, new KrakenResultFastqMergeListener() {
					private long lastLine = -1;

					@Override
					public void newTaxIdForRead(long lineCount, byte[] readDescriptor, byte[] read, byte[] readProbs,
							String krakenTaxid, int bps, String kmerTaxid, int hitLength, byte[] output) {
						if (taxIds == null || taxIds.contains(kmerTaxid)) {
							kmerCountTrie.add(kmerTaxid, hitLength);
							System.out.println(ByteArrayToString.toString(output));
						}
						if (lineCount != lastLine) {
							lastLine = lineCount;
							if (taxIds == null || taxIds.contains(krakenTaxid)) {
								classCountTrie.inc(krakenTaxid);
							}
						}
					}
				});
			}
		};
		if (getLogger().isInfoEnabled()) {
			String execLine = krakenExecutor.genExecLine(getProject().getKrakenDB(), getProject().getFastqFile());
			getLogger().info("Run kraken with " + execLine);
		}
		try {
			krakenExecutor.execute(getProject().getKrakenDB(), getProject().getFastqFile(), null);

			PrintStream out = new PrintStream(new FileOutputStream(file));

			Map<String, Long> map = new HashMap<String, Long>();
			out.println("taxid;kmers");
			kmerCountTrie.collect(map);
			CountingDigitTrie.print(map, out);

			map.clear();
			out.println("taxid;classifications");
			classCountTrie.collect(map);
			CountingDigitTrie.print(map, out);

			out.close();
		} catch (InterruptedException | IOException e) {
			throw new RuntimeException(e);
		}
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Finished kraken");
		}
	}
}
