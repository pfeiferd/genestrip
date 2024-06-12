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
package org.metagene.genestrip.goals.genbank;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.genbank.AssemblySummaryReader;
import org.metagene.genestrip.genbank.AssemblySummaryReader.FTPEntryQuality;
import org.metagene.genestrip.genbank.AssemblySummaryReader.FTPEntryWithQuality;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class FastaFilesFromGenbankGoal extends ObjectGoal<Map<TaxIdNode, List<FTPEntryWithQuality>>, GSProject> {
	private static final Comparator<FTPEntryWithQuality> COMPARATOR = new Comparator<FTPEntryWithQuality>() {
		@Override
		public int compare(FTPEntryWithQuality o1, FTPEntryWithQuality o2) {
			// Lower quality first in resulting sort order...
			return o2.getQuality().ordinal() - o1.getQuality().ordinal();
		}
	};
	
	private final FileGoal<GSProject> assemblyGoal;
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxidsFromGenbankGoal;
	private final List<FTPEntryQuality> fastaQualities;
	private final int maxFromGenbank;

	@SuppressWarnings("unchecked")
	public FastaFilesFromGenbankGoal(GSProject project, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			FileGoal<GSProject> assemblyGoal, ObjectGoal<Set<TaxIdNode>, GSProject> taxidsFromGenbankGoal) {
		super(project, GSGoalKey.FASTAGSENBANK, taxTreeGoal, taxidsFromGenbankGoal, assemblyGoal);
		this.assemblyGoal = assemblyGoal;
		this.taxTreeGoal = taxTreeGoal;
		this.taxidsFromGenbankGoal = taxidsFromGenbankGoal;
		this.fastaQualities = (List<FTPEntryQuality>) configValue(GSConfigKey.FASTA_QUALITIES);
		this.maxFromGenbank = intConfigValue(GSConfigKey.MAX_FROM_GENBANK);
	}

	@Override
	public void makeThis() {
		try {
			if (getLogger().isInfoEnabled()) {
				getLogger().info(
						"Tax ids used for additional fasta downloads from genbank: " + taxidsFromGenbankGoal.get());
			}
			Map<TaxIdNode, List<FTPEntryWithQuality>> res = new HashMap<TaxIdNode, List<FTPEntryWithQuality>>();
			if (!taxidsFromGenbankGoal.get().isEmpty()) {
				// Optimization: Download the assembly summary file only if necessary.
				assemblyGoal.make();

				AssemblySummaryReader assemblySummaryReader = new AssemblySummaryReader(
						getProject().getCommon().getGenbankDir(), true, taxTreeGoal.get());
				int[] nEntriesTotal = new int[1];
				Map<TaxIdNode, List<FTPEntryWithQuality>> entries = assemblySummaryReader.getRelevantEntries(
						taxidsFromGenbankGoal.get(), fastaQualities, nEntriesTotal);
				if (getLogger().isInfoEnabled()) {
					int sum = 0;
					for (TaxIdNode node : entries.keySet()) {
						sum += entries.get(node).size();
					}
					getLogger().info("Potentially relevant entries: " + sum);
				}
				for (TaxIdNode node : entries.keySet()) {
					for (FTPEntryWithQuality entry : entries.get(node)) {
						List<FTPEntryWithQuality> list = res.get(node);
						if (entry != null && isMatchingEntryForNode(entry, node)) {
							if (list == null) {
								list = new ArrayList<FTPEntryWithQuality>();
								res.put(node, list);
							}
							if (!list.contains(entry)) {
								list.add(entry);
								updateEntriesForNode(list, node);
							}
						}
					}					
				}
				// Reduce the number of entries to maxFromGenbank and keep the best ones in terms of quality.
				for (TaxIdNode node : res.keySet()) {
					List<FTPEntryWithQuality> list = res.get(node);
					if (maxFromGenbank > 0 && list.size() > maxFromGenbank) {
						Collections.sort(list, COMPARATOR);
						while (list.size() > maxFromGenbank) {
							list.remove(0);
						}
					}
				}

				if (getLogger().isInfoEnabled()) {
					int sum = 0;
					for (TaxIdNode node : res.keySet()) {
						sum += res.get(node).size();
					}
					getLogger().info("Entries selected for download: " + sum);
					getLogger().info("Total number of entries in assembly summary file: " + nEntriesTotal[0]);
				}
			}
			set(res);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected boolean isMatchingEntryForNode(FTPEntryWithQuality entry, TaxIdNode node) {
		return fastaQualities.contains(entry.getQuality());
	}

	// Just for potential override.
	protected void updateEntriesForNode(List<FTPEntryWithQuality> currentEntries, TaxIdNode node) {
	}
}
