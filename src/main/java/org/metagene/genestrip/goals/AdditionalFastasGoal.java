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
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.genbank.AssemblySummaryReader.FTPEntryWithQuality;
import org.metagene.genestrip.goals.genbank.FastaFilesGenbankDownloadGoal;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class AdditionalFastasGoal extends ObjectGoal<Map<File, TaxIdNode>, GSProject> {
	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#')
			.setDelimiter(' ').setRecordSeparator('\n').build();

	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final ObjectGoal<Map<TaxIdNode, List<FTPEntryWithQuality>>, GSProject> fromGenbank;
	private final FastaFilesGenbankDownloadGoal fastaFilesGenbankDownloadGoal;

	@SafeVarargs
	public AdditionalFastasGoal(GSProject project, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			ObjectGoal<Map<TaxIdNode, List<FTPEntryWithQuality>>, GSProject> fromGenbank,
			FastaFilesGenbankDownloadGoal fastaFilesGenbankDownloadGoal, Goal<GSProject>... dependencies) {
		super(project, GSGoalKey.ADD_FASTAS, append(dependencies, taxTreeGoal, fromGenbank, fastaFilesGenbankDownloadGoal));
		this.taxTreeGoal = taxTreeGoal;
		this.fromGenbank = fromGenbank;
		this.fastaFilesGenbankDownloadGoal = fastaFilesGenbankDownloadGoal;
	}

	@Override
	protected void doMakeThis() {
		Map<File, TaxIdNode> res = new HashMap<File, TaxTree.TaxIdNode>();
		File additonalEntryFile = getProject().getAdditionalFile();
		if (additonalEntryFile.exists()) {
			try (CSVParser parser = FORMAT
					.parse(new InputStreamReader(StreamProvider.getInputStreamForFile(additonalEntryFile)))) {
				TaxTree taxTree = taxTreeGoal.get();
				for (CSVRecord record : parser) {
					String taxid = record.get(0);
					String fastaFilePath = record.get(1);
					TaxIdNode node = taxTree.getNodeByTaxId(taxid);
					if (node != null) {
						File file = getProject().fastaFileFromPath(fastaFilePath);
						if (file != null) {
							res.put(file, node);
						} else if (getLogger().isWarnEnabled()) {
							getLogger().warn("Ignoring missing file " + fastaFilePath + ".");
						}
					} else {
						if (getLogger().isWarnEnabled()) {
							getLogger()
									.warn("Unknown taxid in additional file (omitting fasta files for it): " + taxid);
						}
					}
				}
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
		for (TaxIdNode node : fromGenbank.get().keySet()) {
			for (FTPEntryWithQuality entry : fromGenbank.get().get(node)) {
				File file = fastaFilesGenbankDownloadGoal.entryToFile(entry);
				res.put(file, node);
			}
		}

		set(res);
	}
}
