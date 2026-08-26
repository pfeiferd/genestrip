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
package org.metagene.genestrip.goals.refseq;

import java.util.List;
import java.util.Set;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSConfigKey.SeqType;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.refseq.AccessionFileProcessor;
import org.metagene.genestrip.refseq.AccessionMap;
import org.metagene.genestrip.refseq.AccessionMapTrieImpl;
import org.metagene.genestrip.refseq.RefSeqCategory;
import org.metagene.genestrip.refseq.GenomeKeyTrie;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;

import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxNodeSelection;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * Goal that builds the {@link AccessionMap}: it scans the RefSeq release catalog and maps every
 * accession range to its taxonomy node, incrementing that node's RefSeq contig count. Depends on
 * the tax tree, the selected categories and the downloaded catalog.
 *
 * @param <P> the project type
 */
public class AccessionMapGoal<P extends GSProject> extends ObjectGoal<AccessionMap, P> implements Goal.LogHeapInfo {
	private final ObjectGoal<TaxTree, P> taxTreeGoal;
	private final RefSeqCatalogDownloadGoal catalogGoal;
	private final ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal;
	private final ObjectGoal<TaxNodeSelection, P> taxNodesGoal;

	/**
	 * Creates the goal, wiring the tax tree, categories and catalog download goals it reads.
	 *
	 * @param project the project
	 * @param categoriesGoal goal providing the selected RefSeq categories
	 * @param taxTreeGoal goal providing the taxonomy tree
	 * @param taxNodesGoal goal providing the requested tax nodes, whose genomes are the ones a limit counts
	 * @param catalogGoal goal providing the downloaded RefSeq catalog
	 * @param deps additional goal dependencies
	 */
	@SafeVarargs
	public AccessionMapGoal(P project, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
			ObjectGoal<TaxTree, P> taxTreeGoal, ObjectGoal<TaxNodeSelection, P> taxNodesGoal,
			RefSeqCatalogDownloadGoal catalogGoal, Goal<P>... deps) {
		super(project, GSGoalKey.ACCMAP,
				Goal.append(deps, taxTreeGoal, taxNodesGoal, categoriesGoal, catalogGoal));
		this.categoriesGoal = categoriesGoal;
		this.taxTreeGoal = taxTreeGoal;
		this.taxNodesGoal = taxNodesGoal;
		this.catalogGoal = catalogGoal;
	}

	@Override
	protected void doMakeThis() {
		AccessionFileProcessor processor = new AccessionFileProcessor(categoriesGoal.get(),
				(SeqType) configValue(GSConfigKey.SEQ_TYPE), (List<GSConfigKey.RefSeqStatus>) configValue(GSConfigKey.RES_SEQ_STATUS)) {
			private TaxTree taxTree = taxTreeGoal.get();
			// The map grows as entries arrive, so the catalog is read once here instead of once to
			// count the entries and once more to fill them in.
			// The selection of genomes maxGenomesPerTaxid allows, made here because here it can be made
			// once and the same way every time: one pass, one thread, a fixed order. Null while the
			// limit is not set, and then nothing of this runs. Declared before the map, which carries it
			// from its construction on, so that the map is never seen without the selection it belongs to.
			private final int maxGenomes = intConfigValue(GSConfigKey.MAX_GENOMES_PER_TAXID);
			private final Rank limitRank = (Rank) configValue(GSConfigKey.MAX_PER_TAXID_RANK);
			private final boolean assemblyOnly = booleanConfigValue(GSConfigKey.ASSEMBLY_ACCESSIONS_ONLY);
			private final GenomeKeyTrie admittedGenomes =
					maxGenomes == Integer.MAX_VALUE ? null : new GenomeKeyTrie();
			private final AccessionMap map = new AccessionMapTrieImpl(admittedGenomes);
			// How many genomes each limit node has taken, kept only while the selection is being made.
			// Keyed by the node rather than by its tax id: the node is in hand, its identity hash is one
			// read, and a digit trie over tax id strings would buy nothing here - the key is not a byte
			// range in a buffer, which is the one thing such a trie is for.
			private final Object2IntOpenHashMap<TaxIdNode> genomesPerNode = new Object2IntOpenHashMap<>();
			private final Set<TaxIdNode> requested = taxNodesGoal.get().getSelected();

			@Override
			public void processCatalog(StreamingResource catalogFile) {
				super.processCatalog(catalogFile);
				if (admittedGenomes != null) {
					// The selection is complete here and read from here on, by every pass and every
					// reader thread; closing it turns a later addition into a failure instead of a race.
					admittedGenomes.freeze();
				}
				map.optimize();
				set(map);
			}

			@Override
			protected void handleEntry(byte[] target, int taxIdEnd, int accessionStart, int accessionEnd) {
				TaxIdNode node = taxTree.getNodeByTaxId(target, 0, taxIdEnd);
				if (node != null) {
					map.put(target, accessionStart, accessionEnd, node);
					node.incRefSeqContigs();
					admit(target, accessionStart, accessionEnd, node);
				}
			}

			/**
			 * Takes the genome of this accession into the selection if its taxon still has room, and
			 * does nothing at all while no limit is set.
			 * <p>
			 * A genome counts once however many accessions it arrives in, which is what the set is for -
			 * the contigs of one assembly need not be contiguous in the catalog for that to hold. Only
			 * the requested taxa are counted: a taxon whose contigs the fill never stores must not take
			 * places from one whose contigs it does. And the limit node is the ancestor at
			 * {@code maxPerTaxidRank}, or the node itself where the lineage has none.
			 */
			private void admit(byte[] target, int accessionStart, int accessionEnd, TaxIdNode node) {
				if (admittedGenomes == null || !countsAsGenome(target, accessionStart)
						|| !(requested.isEmpty() || requested.contains(node))
						|| admittedGenomes.isAdmitted(target, accessionStart, accessionEnd)) {
					return;
				}
				TaxIdNode limitNode = limitNodeOf(node);
				if (genomesPerNode.getInt(limitNode) < maxGenomes) {
					admittedGenomes.admit(target, accessionStart, accessionEnd);
					genomesPerNode.addTo(limitNode, 1);
				}
			}

			/**
			 * Returns whether this accession is one of a genome, applying the same filter the map's
			 * lookups apply - a limit counted over accessions the fill never resolves would be counted
			 * against the wrong denominator.
			 */
			private boolean countsAsGenome(byte[] target, int accessionStart) {
				return assemblyOnly ? isAssemblyAccession(target, accessionStart)
						: isGenomicAccession(target, accessionStart);
			}

			/** The node a limit is counted at: the ancestor at {@code maxPerTaxidRank}, or the node. */
			private TaxIdNode limitNodeOf(TaxIdNode node) {
				if (limitRank != null) {
					for (TaxIdNode n = node; n != null; n = n.getParent()) {
						if (limitRank.equals(n.getRank())) {
							return n;
						}
					}
				}
				return node;
			}

			@Override
			protected boolean isProgressBar() {
				return booleanConfigValue(GSConfigKey.PROGRESS_BAR);
			}

			protected String getProgressBarTaskName() {
				return getKey().getName();
			}
		};
		processor.processCatalog(new StreamingFileResource(catalogGoal.getCatalogFile()));
	}
}
