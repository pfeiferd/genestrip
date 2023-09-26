package org.metagene.genestrip.goals;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class AdditionalFastasGoal extends ObjectGoal<Map<File, TaxIdNode>, GSProject> {
	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#')
			.setDelimiter('\t').setRecordSeparator('\n').build();

	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;

	@SafeVarargs
	public AdditionalFastasGoal(GSProject project, String name, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			Goal<GSProject>... dependencies) {
		super(project, name, append(dependencies, taxTreeGoal));
		this.taxTreeGoal = taxTreeGoal;
	}

	@Override
	public void makeThis() {
		try {
			Map<File, TaxIdNode> res = new HashMap<File, TaxTree.TaxIdNode>();

			File additonalEntryFile = getProject().getAddtionalFile();
			if (additonalEntryFile.exists()) {
				TaxTree taxTree = taxTreeGoal.get();

				Reader in = new InputStreamReader(StreamProvider.getInputStreamForFile(additonalEntryFile));

				Iterable<CSVRecord> records = FORMAT.parse(in);
				for (CSVRecord record : records) {
					String taxid = record.get(0);
					String fileName = record.get(1);

					TaxIdNode node = taxTree.getNodeByTaxId(taxid);
					if (node != null) {
						File file = new File(getProject().getFastasDir(), fileName);
						if (file.exists()) {
							if (getLogger().isInfoEnabled()) {
								getLogger()
										.warn("Adding additional fasta file " + file + " for taxid " + node.getTaxId());
							}
							res.put(file, node);
						} else {
							if (getLogger().isWarnEnabled()) {
								getLogger().warn("Missing additional fasta file " + file);
							}
						}
					}
				}
			}
			set(res);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
