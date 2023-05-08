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
import java.io.OutputStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.fastqgen.KMerFastqGenerator;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryQuality;
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryWithQuality;
import org.metagene.genestrip.tax.TaxTree.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;

public class KMerFastqGoal extends FileListGoal<GSProject> {
	private final Map<TaxIdNode, Set<File>> taxToFastas;

	@SafeVarargs
	public KMerFastqGoal(GSProject project, String name,
			ObjectGoal<Map<TaxIdNode, List<FTPEntryWithQuality>>, GSProject> fastaFilesGoal, Goal<GSProject>... deps) {
		super(project, name, project.getOutputFile(name, FileType.FASTQ), ArraysUtil.append(deps, fastaFilesGoal));

		Rank maxRank = project.getConfig().getMaxRankForFilters();
		taxToFastas = new HashMap<TaxIdNode, Set<File>>();
		for (TaxIdNode key : fastaFilesGoal.get().keySet()) {
			List<FTPEntryWithQuality> entries = fastaFilesGoal.get().get(key);
			while (key.getRank() != null && key.getRank().isBelow(maxRank)) {
				key = key.getParent();
			}
			if (key.getRank() == null) {
				if (getLogger().isWarnEnabled()) {
					getLogger().warn("Missing rank for taxid: " + key.getName() + " - omitted");
				}				
				continue;
			}
			Set<File> files = taxToFastas.get(key);
			boolean isNew = false;
			if (files == null) {
				isNew = true;
				files = new HashSet<File>();
			}
			FTPEntryQuality minQuality = getProject().getConfig().getFastaQuality();
			for (FTPEntryWithQuality entry : entries) {
				if (minQuality == null || !entry.getQuality().below(minQuality)) {
					files.add(new File(getProject().getFastasDir(), entry.getFileName()));
				}
			}
			if (isNew && !files.isEmpty()) {
				taxToFastas.put(key, files);				
			}
		}
	}

	@Override
	protected void makeFile(File fastq) {
		try {
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Generating fastq file " + fastq);
			}
			OutputStream output = StreamProvider.getOutputStreamForFile(fastq);

			KMerFastqGenerator generator = new KMerFastqGenerator(getName(), getProject().getkMserSize(),
					getProject().getConfig().getKmerFastInitialBloomSize(),
					getProject().getConfig().getKmerFastBloomFpp(), getProject().getConfig().getMaxReadSizeBytes(),
					output, getProject().isIgnoreMissingFiles());

			/*
			 * We run the kmer generation "bit by bit" in chunks by common max rank (usually
			 * genus or species). If we ran the kmers altogether, then the bloom filter in
			 * the background might (inside KMerFastqGenerator) would consume way too much
			 * memory for larger projects (and would fail). We need the bloom filter inside
			 * KMerFastqGenerator in order to ensure uniqueness (no double entries) for
			 * k-mers for taxons on the max rank level or below. We need this to make sure
			 * that the finally generated filter has the right size (in terms of entries).
			 * It will only contain k-mers on max rank level or below anyways so the
			 * uniqueness on this level is sufficient here.
			 */
			int max = taxToFastas.size();
			int counter = 0;
			for (TaxIdNode taxId : taxToFastas.keySet()) {
				counter++;
				if (getLogger().isInfoEnabled()) {
					getLogger().info("Processing taxid (" + counter + "/" + max + "): " + taxId.getName());
				}
				long addedKmers = generator.add(taxToFastas.get(taxId));
				if (getLogger().isInfoEnabled()) {
					getLogger().info("Entered K-mers: " + addedKmers);
				}
			}

			output.close();
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Generated fastq file " + fastq);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
}
