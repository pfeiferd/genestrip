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
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.Rank;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;

public class KMerFastqGoal extends FileListGoal<GSProject> {
	private final ObjectGoal<Map<TaxIdNode, List<FTPEntryWithQuality>>, GSProject> fastaFilesGoal;
	private final KMerFastqGenerator generator;
	private final File tempFile;
	
	private Map<File, TaxIdNode> fileToCollectNode;
	private Map<TaxIdNode, Map<TaxIdNode,File>> taxToFastas;

	@SafeVarargs
	public KMerFastqGoal(GSProject project, String name,
			ObjectGoal<Map<TaxIdNode, List<FTPEntryWithQuality>>, GSProject> fastaFilesGoal, Goal<GSProject>... deps) {
		super(project, name, (List<File>) null, ArraysUtil.append(deps, fastaFilesGoal));
		this.fastaFilesGoal = fastaFilesGoal;
		generator = new KMerFastqGenerator(getName(), getProject().getKMserSize(),
				getProject().getConfig().getKMerFastBloomFpp(), getProject().getConfig().getMaxReadSizeBytes(),
				getProject().isIgnoreMissingFiles());

		tempFile = project.getOutputFile(name + "_temp", FileType.FASTQ);
	}
	
	@Override
	protected void provideFiles() {
		Rank maxRank = getProject().getConfig().getMaxRankForFilters();
		taxToFastas = new HashMap<TaxIdNode, Map<TaxIdNode,File>>();
		for (TaxIdNode key : fastaFilesGoal.get().keySet()) {
			TaxIdNode fastaFileNode = key;
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
			Map<TaxIdNode, File> files = taxToFastas.get(key);
			boolean isNew = false;
			if (files == null) {
				isNew = true;
				files = new HashMap<TaxIdNode,File>();
			}
			FTPEntryQuality minQuality = getProject().getConfig().getFastaQuality();
			for (FTPEntryWithQuality entry : entries) {
				if (minQuality == null || !entry.getQuality().below(minQuality)) {
					files.put(fastaFileNode, new File(getProject().getFastasDir(), entry.getFileName()));
				}
			}
			if (isNew && !files.isEmpty()) {
				taxToFastas.put(key, files);
			}
		}
		Set<TaxIdNode> aboveMaxRank = new HashSet<TaxIdNode>();
		for (TaxIdNode taxId : taxToFastas.keySet()) {
			aboveMaxRank.add(taxId.getParent());
		}
		fileToCollectNode = new HashMap<File, TaxTree.TaxIdNode>();
		for (TaxIdNode collectNode : aboveMaxRank) {
			File file = getProject().getOutputFile(getName() + "_" + collectNode.getTaxId(), FileType.FASTQ);
			fileToCollectNode.put(file, collectNode);
			addFile(file);
		}
	}
	
	@Override
	public List<File> getFiles() {
		return super.getFiles();
	}

	@Override
	protected void makeFile(File fastq) {
		try {
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Generating fastq file " + fastq);
			}

			TaxIdNode collectNode = fileToCollectNode.get(fastq);

			if (tempFile.exists()) {
				tempFile.delete();
			}

			OutputStream output = StreamProvider.getOutputStreamForFile(tempFile);
			generator.setOutputStream(output);

			int counter = 0;
			List<TaxIdNode> subnodes = collectNode.getSubNodes();
			int max = subnodes.size();

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
			for (TaxIdNode child : subnodes) {
				counter++;
				Map<TaxIdNode, File> fastas = taxToFastas.get(child);
				if (fastas != null) {
					if (getLogger().isInfoEnabled()) {
						getLogger().info("Processing taxid (" + counter + "/" + max + "): " + child.getName());
					}
					long addedKmers = generator.add(fastas);
					if (getLogger().isInfoEnabled()) {
						getLogger().info("Entered K-mers: " + addedKmers);
					}
				}
			}

			output.close();

			tempFile.renameTo(fastq);

			if (getLogger().isInfoEnabled()) {
				getLogger().info("Generated fastq file " + fastq);
			}
		} catch (

		IOException e) {
			throw new RuntimeException(e);
		}
	}
}
