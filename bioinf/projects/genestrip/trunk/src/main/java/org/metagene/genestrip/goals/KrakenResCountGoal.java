package org.metagene.genestrip.goals;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Map;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.kraken.KrakenExecutor;
import org.metagene.genestrip.kraken.KrakenResultParser;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.util.CountingDigitTrie;

public class KrakenResCountGoal extends FileListGoal<GSProject> {
	@SafeVarargs
	public KrakenResCountGoal(GSProject project, Goal<GSProject>... dependencies) {
		super(project, "krakenrescount", project.getTaxCountsFile(), dependencies);
	}
	
	@Override
	protected void makeFile(File file) throws IOException {
		@SuppressWarnings("unchecked")
		Map<String, Long>[] res = new Map[1];
		
		KrakenExecutor krakenExecutor = new KrakenExecutor(getProject().getConfig().getKrakenBinFolder(),
				getProject().getConfig().getKrakenExecExpr()) {
			@Override
			protected void handleOutputStream(InputStream stream, File outFile) throws IOException {
				KrakenResultParser parser = new KrakenResultParser();
				
				res[0] = parser.process(new BufferedInputStream(stream), null);
			}
		};
		if (getLogger().isInfoEnabled()) {
			String execLine = krakenExecutor.genExecLine(getProject().getKrakenDB(), getProject().getFastqFile());
			getLogger().info("Run kraken with " + execLine);
		}
		try {
			krakenExecutor.execute(getProject().getKrakenDB(), getProject().getFastqFile(), null);
			FileOutputStream out = new FileOutputStream(file);
			CountingDigitTrie.print(res[0], out);
			out.close();
		} catch (InterruptedException | IOException e) {
			throw new RuntimeException(e);
		}
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Finished kraken");
		}
	}
}
