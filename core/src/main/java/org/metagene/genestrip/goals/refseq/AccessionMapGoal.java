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

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import org.metagene.genestrip.genbank.AssemblySummaryReader;
import org.metagene.genestrip.refseq.ReferenceGenomes;
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
import org.metagene.genestrip.refseq.AccessionTrie;
import org.metagene.genestrip.refseq.AssemblyInfo;
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
	private final AssemblyMetadataGoal<P> assemblyMetaGoal;
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
	 * @param assemblyMetaGoal the goal providing, for every sequence of an admitted assembly, the key
	 *            of that assembly; consulted only where {@code refseq.genomesOnly} is not off.
	 * @param deps additional goal dependencies
	 */
	@SafeVarargs
	public AccessionMapGoal(P project, ObjectGoal<Set<RefSeqCategory>, P> categoriesGoal,
			ObjectGoal<TaxTree, P> taxTreeGoal, ObjectGoal<TaxNodeSelection, P> taxNodesGoal,
			RefSeqCatalogDownloadGoal catalogGoal, AssemblyMetadataGoal<P> assemblyMetaGoal, Goal<P>... deps) {
		super(project, GSGoalKey.ACCMAP,
				Goal.append(deps, taxTreeGoal, taxNodesGoal, categoriesGoal, catalogGoal));
		this.categoriesGoal = categoriesGoal;
		this.taxTreeGoal = taxTreeGoal;
		this.taxNodesGoal = taxNodesGoal;
		this.catalogGoal = catalogGoal;
		this.assemblyMetaGoal = assemblyMetaGoal;
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
			private final GSConfigKey.GenomesOnly genomesOnly =
					(GSConfigKey.GenomesOnly) configValue(GSConfigKey.GENOMES_ONLY);
			// `preferRefFirst' restricts nothing -- drafts enter as they do while the filter is off --
			// and only decides which genomes a limit lets in. Every other setting but off restricts to
			// whole assemblies; which levels count as whole is decided where the index is built, so
			// that whatever it holds is something to keep.
			private final boolean preferReference = genomesOnly == GSConfigKey.GenomesOnly.PREFER_REF_FIRST;
			private final boolean completeOnly =
					genomesOnly != GSConfigKey.GenomesOnly.OFF && !preferReference;
			// The assembly metadata answers two different questions and is read for both: which
			// sequences belong to a whole assembly, and which assembly a species is represented by.
			private final boolean withAssemblyInfo = completeOnly || preferReference;
			// Asked for only where the filter is on, which is what keeps the metadata goal from being
			// made - and so the assembly summary from being read and the catalog from being walked a
			// second time - in a build that does not want it.
			private final AccessionTrie<AssemblyInfo> completeGenomes = withAssemblyInfo ? assemblyMetaGoal.get() : null;
			// The marked reference assemblies, read once here: which taxa have one, so that a place can
			// be kept, and the genome keys of the marked drafts, which the length index cannot name.
			private final ReferenceGenomes references = preferReference ? readReferences() : null;
			// The limit nodes whose marked assembly has arrived. Until it has, one of their places is
			// held back; afterwards they fill like any other.
			private final Set<TaxIdNode> referenceAdmitted = preferReference ? new HashSet<>() : null;
			private final GenomeKeyTrie admittedGenomes =
					maxGenomes == Integer.MAX_VALUE ? null : new GenomeKeyTrie();
			private final AccessionMap map = new AccessionMapTrieImpl(admittedGenomes);
			// How many genomes each limit node has taken, kept only while the selection is being made.
			// Keyed by the node rather than by its tax id: the node is in hand, its identity hash is one
			// read, and a digit trie over tax id strings would buy nothing here - the key is not a byte
			// range in a buffer, which is the one thing such a trie is for. Both of these exist only
			// where a selection is being made; without a limit nothing of it is built or asked for.
			private final Object2IntOpenHashMap<TaxIdNode> genomesPerNode =
					admittedGenomes == null ? null : new Object2IntOpenHashMap<>();
			private final Set<TaxIdNode> requested =
					admittedGenomes == null ? null : taxNodesGoal.get().getSelected();

			/**
			 * Reads the marked reference assemblies off the assembly summary, which is where the
			 * metadata goal takes its own index from. Read here rather than carried by that goal: it
			 * is one scan of a file already on disk, and the goal's result is an accession trie whose
			 * shape has nothing to say about which assembly a species is represented by.
			 */
			private ReferenceGenomes readReferences() {
				File summary = new File(getProject().getCommon().getRefSeqDir(),
						AssemblySummaryReader.ASSEMLY_SUM_REFSEQ);
				try {
					ReferenceGenomes result = new ReferenceGenomes(summary);
					if (getLogger().isInfoEnabled()) {
						getLogger().info("Reference genomes: " + result.taxonCount() + " taxa marked, "
								+ result.draftCount() + " of them by a draft assembly.");
					}
					return result;
				} catch (IOException e) {
					throw new RuntimeException("Cannot read " + summary + ", which "
							+ GSConfigKey.GENOMES_ONLY.getName() + "="
							+ GSConfigKey.GenomesOnly.PREFER_REF_FIRST.getName() + " needs.", e);
				}
			}

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
				// The sequences of a draft are kept out of the map rather than filtered at lookup, so that
				// the fill, the genome cap and the reported contig counts all speak about the same ones.
				// What counts as finished is not guessed from the accession: the metadata goal has matched
				// the catalog's sequences against NCBI's assembly summary and says which belong to an
				// assembly that summary calls complete. RNA accessions pass unchanged, as they do for
				// refseq.assemblyAccessionsOnly: the filter is about assemblies and must not be able to
				// empty a database of transcripts.
				if (completeOnly && isGenomicAccession(target, accessionStart)
						&& assemblyOf(target, accessionStart, accessionEnd) == null) {
					return;
				}
				TaxIdNode node = taxTree.getNodeByTaxId(target, 0, taxIdEnd);
				if (node != null) {
					map.put(target, accessionStart, accessionEnd, node);
					node.incRefSeqContigs();
					if (admittedGenomes != null) {
						admit(target, accessionStart, accessionEnd, node);
					}
				}
			}

			/**
			 * Takes the genome of this accession into the selection if its taxon still has room.
			 * <p>
			 * A genome counts once however many accessions it arrives in, which is what the set is for -
			 * the contigs of one assembly need not be contiguous in the catalog for that to hold. Only
			 * the requested taxa are counted: a taxon whose contigs the fill never stores must not take
			 * places from one whose contigs it does. And the limit node is the ancestor at
			 * {@code maxPerTaxidRank}, or the node itself where the lineage has none.
			 * <p>
			 * Called only where a selection is being made; the caller checks.
			 */
			private void admit(byte[] target, int accessionStart, int accessionEnd, TaxIdNode node) {
				if (!countsAsGenome(target, accessionStart)
						|| !(requested.isEmpty() || requested.contains(node))
						|| isAlreadyAdmitted(target, accessionStart, accessionEnd)) {
					return;
				}
				TaxIdNode limitNode = limitNodeOf(node);
				// Counted by assembly where one is known. The genome key groups the contigs of a
				// shotgun project but leaves the chromosome and each plasmid of a finished genome
				// standing for themselves, so a limit counted on it would admit a quarter of what it
				// was asked for wherever finished genomes carry plasmids.
				AssemblyInfo info = withAssemblyInfo ? assemblyOf(target, accessionStart, accessionEnd) : null;
				boolean marked = preferReference && isReference(target, accessionStart, accessionEnd, info);
				if (genomesPerNode.getInt(limitNode) < roomFor(limitNode, marked)) {
					byte[] key = info == null ? null : info.getKey();
					if (key != null) {
						admittedGenomes.admitKey(key, 0, key.length);
					} else {
						admittedGenomes.admit(target, accessionStart, accessionEnd);
					}
					genomesPerNode.addTo(limitNode, 1);
					if (marked) {
						referenceAdmitted.add(limitNode);
					}
				}
			}

			/**
			 * How many genomes this limit node may hold before the genome in hand is turned away. It is
			 * the limit itself, except that a taxon whose reference assembly is still to come keeps one
			 * place back for it: without that the limit would be full by the time the marked assembly
			 * arrives, and which genomes got in would again be a matter of the order the catalog states
			 * them in. The place is held for the marked assembly alone, so a taxon whose marked
			 * assembly the release does not carry ends one genome short -- the price of holding the
			 * limit exactly rather than exceeding it by one.
			 */
			private int roomFor(TaxIdNode limitNode, boolean marked) {
				if (!preferReference) {
					return maxGenomes;
				}
				return ReferenceGenomes.roomFor(maxGenomes, marked, referenceAdmitted.contains(limitNode),
						references.hasReference(limitNode.getTaxId()));
			}

			/**
			 * Whether this accession's genome is the assembly its taxon is represented by. A finished
			 * one is known from the length index, which has matched it to its sequences; a draft is not
			 * in that index and is recognised by the WGS prefix its contigs carry.
			 */
			private boolean isReference(byte[] target, int accessionStart, int accessionEnd, AssemblyInfo info) {
				return (info != null && info.isReference())
						|| references.isReferenceDraft(target, accessionStart, accessionEnd);
			}

			/**
			 * Returns whether this accession is one of a genome, applying the same filter the map's
			 * lookups apply - a limit counted over accessions the fill never resolves would be counted
			 * against the wrong denominator.
			 */
			/**
			 * Whether this accession's genome - its assembly where one is known - is already in. Asked
			 * with the same key the admission uses, or a genome counted by assembly there and by
			 * accession here would be taken in once per replicon.
			 */
			private boolean isAlreadyAdmitted(byte[] target, int accessionStart, int accessionEnd) {
				AssemblyInfo info = withAssemblyInfo ? assemblyOf(target, accessionStart, accessionEnd) : null;
				byte[] key = info == null ? null : info.getKey();
				return key != null ? admittedGenomes.isAdmittedKey(key, 0, key.length)
						: admittedGenomes.isAdmitted(target, accessionStart, accessionEnd);
			}

			/**
			 * Returns the key of the assembly this accession belongs to, or null where it belongs to none
			 * the summary calls complete. Null is also the answer while the filter is off, which no caller
			 * asks for.
			 */
			private AssemblyInfo assemblyOf(byte[] target, int accessionStart, int accessionEnd) {
				int keyEnd = accessionStart
						+ GenomeKeyTrie.genomeKeyLength(target, accessionStart, accessionEnd);
				return completeGenomes.get(target, accessionStart, keyEnd);
			}

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
