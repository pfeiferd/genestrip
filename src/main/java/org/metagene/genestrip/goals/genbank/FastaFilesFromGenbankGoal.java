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
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.genbank.AssemblySummaryReader;
import org.metagene.genestrip.genbank.AssemblySummaryReader.AssemblyQuality;
import org.metagene.genestrip.genbank.AssemblySummaryReader.AssemblyEntry;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

public class FastaFilesFromGenbankGoal extends ObjectGoal<Map<TaxIdNode, List<AssemblyEntry>>, GSProject> {
	private static final Comparator<AssemblyEntry> COMPARATOR = new Comparator<AssemblyEntry>() {
		@Override
		public int compare(AssemblyEntry o1, AssemblyEntry o2) {
			// Lower quality first in resulting sort order...
			return o2.getQuality().ordinal() - o1.getQuality().ordinal();
		}
	};
	
	private final FileGoal<GSProject> assemblyGoal;
	private final ObjectGoal<TaxTree, GSProject> taxTreeGoal;
	private final ObjectGoal<Set<TaxIdNode>, GSProject> taxidsFromGenbankGoal;
	private final List<AssemblyQuality> fastaQualities;
	private final int maxFromGenbank;
	private final boolean refGenOnly;

	@SuppressWarnings("unchecked")
	public FastaFilesFromGenbankGoal(GSProject project, ObjectGoal<TaxTree, GSProject> taxTreeGoal,
			FileGoal<GSProject> assemblyGoal, ObjectGoal<Set<TaxIdNode>, GSProject> taxidsFromGenbankGoal) {
		super(project, GSGoalKey.FASTAGSENBANK, taxTreeGoal, taxidsFromGenbankGoal, assemblyGoal);
		this.assemblyGoal = assemblyGoal;
		this.taxTreeGoal = taxTreeGoal;
		this.taxidsFromGenbankGoal = taxidsFromGenbankGoal;
		this.fastaQualities = (List<AssemblyQuality>) configValue(GSConfigKey.FASTA_QUALITIES);
		this.maxFromGenbank = intConfigValue(GSConfigKey.MAX_FROM_GENBANK);
		this.refGenOnly = booleanConfigValue(GSConfigKey.REF_GEN_ONLY);
	}

	@Override
	protected void doMakeThis() {
		if (taxidsFromGenbankGoal.get().isEmpty()) {
			set(Collections.emptyMap());
		}
		else {
			try {
				if (getLogger().isDebugEnabled()) {
					getLogger().debug(
							"Tax ids used for additional fasta downloads from genbank: " + taxidsFromGenbankGoal.get());
				}
				// Optimization: Download the assembly summary file only if necessary.
				assemblyGoal.make();

				AssemblySummaryReader assemblySummaryReader = new AssemblySummaryReader(
						getProject().getCommon().getGenbankDir(), true, taxTreeGoal.get());
				int[] nEntriesTotal = new int[1];
				Map<TaxIdNode, List<AssemblyEntry>> entries = assemblySummaryReader.getRelevantEntries(
						taxidsFromGenbankGoal.get(), fastaQualities, refGenOnly, false, nEntriesTotal);
				if (getLogger().isInfoEnabled()) {
					int sum = 0;
					for (TaxIdNode node : entries.keySet()) {
						sum += entries.get(node).size();
					}
					if (getLogger().isInfoEnabled()) {
						getLogger().info("Potentially relevant entries: " + sum);
					}
				}

				// Reduce the number of entries to maxFromGenbank and keep the best ones in terms of quality.
				for (TaxIdNode node : entries.keySet()) {
					List<AssemblyEntry> list = entries.get(node);
					if (maxFromGenbank > 0 && list.size() > maxFromGenbank) {
						Collections.sort(list, COMPARATOR);
						while (list.size() > maxFromGenbank) {
							list.remove(0);
						}
					}
				}

				if (getLogger().isInfoEnabled()) {
					int sum = 0;
					for (TaxIdNode node : entries.keySet()) {
						sum += entries.get(node).size();
					}
					getLogger().info("Entries selected for download: " + sum);
					getLogger().info("Total number of entries in assembly summary file: " + nEntriesTotal[0]);
				}
				set(entries);
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
	}
}
