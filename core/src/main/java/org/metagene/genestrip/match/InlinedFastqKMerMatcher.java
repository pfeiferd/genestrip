/*
 *
 * "Commons Clause" License Condition v1.0
 *
 * The Software is provided to you by the Licensor under the License,
 * as defined below, subject to the following condition.
 *
 * Without limiting other conditions in the License, the grant of rights under the License
 * will not include, and the License does not grant to you, the right to Sell the Software.
 *
 * For purposes of the foregoing, "Sell" means practicing any or all of the rights granted
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
package org.metagene.genestrip.match;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.util.CGAT;

// A more inlined implementation of the FastqKMerMatcher composed by Claude.
// It does not bring much improvement and therefore, it just "sits here" for unknown future purposed.
// I set it to be deprecated...
@Deprecated
public class InlinedFastqKMerMatcher extends FastqKMerMatcher {

	public InlinedFastqKMerMatcher(KMerSortedArray<SmallTaxIdNode> kmerStore, int initialReadSize,
			int maxQueueSize, ExecutionContext bundle, boolean withProbs, int maxKmerResCounts,
			SmallTaxTree taxTree, int maxPaths, double maxReadTaxErrorCount,
			double maxReadClassErrorCount, boolean writeAll, int threshold, String dbMD5) {
		super(kmerStore, initialReadSize, maxQueueSize, bundle, withProbs, maxKmerResCounts,
				taxTree, maxPaths, maxReadTaxErrorCount, maxReadClassErrorCount, writeAll, threshold, dbMD5);
	}

	@Override
	protected boolean matchRead(final MatcherReadEntry entry, final int index) {
		// Cache all hot fields as local finals to reduce field-dereference overhead
		// in the inner loop. The JIT promotes locals to registers more readily than
		// instance fields.
		final byte[] read = entry.read;
		final int[] badPos = entry.badPos;
		final long[] indexPos = entry.indexPos;
		final int[] jumpTable = CGAT.CGAT_JUMP_TABLE;
		final int[] revJumpTable = CGAT.CGAT_REVERSE_JUMP_TABLE;
		final long straightMask = CGAT.SHIFT_FILTERS_STRAIGHT[k];
		final int reverseShift = (int) CGAT.SHIFT_FILTERS_REVERSE[k];

		boolean found = false;
		int prints = 0;
		int readTaxErrorCount = taxTree == null ? -1 : 0;

		SmallTaxIdNode taxIdNode;
		final int max = entry.readSize - k + 1;
		SmallTaxIdNode lastTaxid = null;
		int contigLen = 0;
		CountsPerTaxid stats = null;

		long kmer = -1;
		long reverseKmer = -1;
		int oldIndex = 0;

		for (int i = 0; i < max; i++) {
			if (kmer == -1) {
				// === Inlined kMerToLongStraight ===
				kmer = 0;
				final int initMax = i + k;
				for (int i1 = i; i1 < initMax; i1++) {
					kmer = (kmer << 2) | (kmer >>> -2);
					final int c = jumpTable[read[i1]];
					if (c == -1) {
						kmer = -1;
						badPos[0] = i1;
						break;
					}
					kmer += c;
				}
				if (kmer == -1) {
					oldIndex = i;
					i = badPos[0];
				} else {
					// === Inlined kMerToLongReverse ===
					// No bad-char check: range [i, i+k) was already verified clean above.
					reverseKmer = 0;
					for (int i1 = i + k - 1; i1 >= i; i1--) {
						reverseKmer = (reverseKmer << 2) | (reverseKmer >>> -2);
						reverseKmer += revJumpTable[read[i1]];
					}
				}
			} else {
				// === Inlined nextKMerStraight + nextKMerReverse ===
				// Read the incoming base once; reuse for both directions.
				final byte newBase = read[i + k - 1];
				final int c = jumpTable[newBase];
				if (c == -1) {
					kmer = -1;
					oldIndex = i;
					i += k - 1;
				} else {
					kmer = ((kmer << 2) & straightMask) | (long) c;
					reverseKmer = (reverseKmer >>> 2) | (((long) revJumpTable[newBase]) << reverseShift);
				}
			}

			// Canonical k-mer: the lexicographically larger of forward and reverse complement.
			taxIdNode = kmer == -1 ? INVALID_NODE
					: kmerStore.getLongInlined(kmer > reverseKmer ? kmer : reverseKmer, indexPos);

			if (readTaxErrorCount != -1) {
				if (taxIdNode == null || taxIdNode == INVALID_NODE) {
					readTaxErrorCount++;
					if (maxReadTaxErrorCount >= 0) {
						if ((maxReadTaxErrorCount >= 1 && readTaxErrorCount > maxReadTaxErrorCount)
								|| (readTaxErrorCount > maxReadTaxErrorCount * max)) {
							readTaxErrorCount = -1;
						}
					}
				} else {
					// === Inlined updateReadTaxid ===
					taxTree.incCount(taxIdNode, index, entry.readNo);
					boolean foundPath = false;
					for (int i1 = 0; i1 < entry.usedPaths; i1++) {
						if (taxTree.isAncestorOf(taxIdNode, entry.readTaxIdNode[i1])) {
							entry.readTaxIdNode[i1] = taxIdNode;
							foundPath = true;
							break;
						} else if (taxTree.isAncestorOf(entry.readTaxIdNode[i1], taxIdNode)) {
							foundPath = true;
							break;
						}
					}
					if (!foundPath && entry.usedPaths < maxPaths) {
						entry.readTaxIdNode[entry.usedPaths++] = taxIdNode;
					}
				}
			}

			if (taxIdNode != lastTaxid) {
				if (contigLen > 0) {
					if (out != null) {
						printKrakenStyleOut(entry, lastTaxid, contigLen, prints++);
					}
					if (stats != null) {
						synchronized (stats) {
							stats.contigs++;
							stats.contigLenSquaredSum += contigLen * contigLen;
							if (contigLen > stats.maxContigLen) {
								stats.maxContigLen = contigLen;
								int j = 1;
								for (; j < entry.readDescriptorSize && entry.readDescriptor[j] != ' '; j++) {
									stats.maxContigDescriptor[j - 1] = entry.readDescriptor[j];
								}
								stats.maxContigDescriptor[j - 1] = 0;
							}
						}
					}
					contigLen = 0;
				}
			}

			if (taxIdNode == INVALID_NODE) {
				contigLen += i >= max ? max - oldIndex : i - oldIndex + 1;
			} else {
				contigLen++;
			}
			lastTaxid = taxIdNode;

			if (taxIdNode != null && taxIdNode != INVALID_NODE) {
				// === Inlined getCountsPerTaxid + stats update ===
				// Use the public field storeIndex directly instead of the getter.
				final int vi = taxIdNode.storeIndex;
				stats = statsIndex[vi];
				if (stats == null) {
					synchronized (statsIndex) {
						if (statsIndex[vi] == null) {
							statsIndex[vi] = new CountsPerTaxid(taxIdNode.getLevel(), taxIdNode.getTaxId(), initialReadSize);
						}
						stats = statsIndex[vi];
					}
				}
				synchronized (stats) {
					stats.kmers++;
					found = true;
					if (readNoPerCPerStat[index][vi] != entry.readNo) {
						stats.reads1KMer++;
						readNoPerCPerStat[index][vi] = entry.readNo;
					}
				}
				if (uniqueCounter != null) {
					uniqueCounter.putInlined(indexPos[0]);
				}
			} else {
				stats = null;
			}
		}

		if (contigLen > 0 && out != null) {
			printKrakenStyleOut(entry, lastTaxid, contigLen, prints);
		}

		if (found) {
			if (contigLen > 0) {
				if (stats != null) {
					synchronized (stats) {
						stats.contigs++;
						stats.contigLenSquaredSum += contigLen * contigLen;
						if (contigLen > stats.maxContigLen) {
							stats.maxContigLen = contigLen;
							int j = 1;
							for (; j < entry.readDescriptorSize && entry.readDescriptor[j] != ' '; j++) {
								stats.maxContigDescriptor[j - 1] = entry.readDescriptor[j];
							}
							stats.maxContigDescriptor[j - 1] = 0;
						}
					}
				}
			}
			if (readTaxErrorCount != -1) {
				int ties = 0;
				for (int i = 0; i < entry.usedPaths; i++) {
					short sum = taxTree.sumCounts(entry.readTaxIdNode[i], index, entry.readNo);
					if (sum > entry.counts[0]) {
						entry.counts[0] = sum;
						entry.readTaxIdNode[0] = entry.readTaxIdNode[i];
						ties = 0;
					} else if (sum == entry.counts[0]) {
						ties++;
						entry.counts[ties] = sum;
						entry.readTaxIdNode[ties] = entry.readTaxIdNode[i];
					}
				}
				if (threshold > 1) {
					for (int i = 0; i <= ties; i++) {
						entry.readTaxIdNode[i] = taxTree.lowestNodeWhereSumAboveThreshold(
								entry.readTaxIdNode[i], index, entry.readNo, threshold);
					}
				}
				SmallTaxIdNode node = entry.readTaxIdNode[0];
				for (int i = 1; i <= ties; i++) {
					node = taxTree.getLowestCommonAncestor(node, entry.readTaxIdNode[i]);
				}
				entry.classNode = node;
				if (node == null) {
					return false;
				}
				// For 'readKmers', count k-mers from entry.readTaxIdNode[0], not just the LCA
				// node. In the tie case the k-mers from one involved path solidify the LCA.
				int readKmers = ties > 0
						? taxTree.sumCounts(entry.readTaxIdNode[0], index, entry.readNo)
						: entry.counts[0];
				int classErrC = max - readKmers;
				if (maxReadClassErrorCount < 0
						|| (maxReadClassErrorCount >= 1 && classErrC <= maxReadClassErrorCount)
						|| (classErrC <= maxReadClassErrorCount * max)) {
					double err = ((double) readTaxErrorCount) / max;
					double classErr = ((double) classErrC) / max;
					entry.classNode = node;
					final int vi = node.storeIndex;
					if (vi >= 0) {
						stats = statsIndex[vi];
						if (stats == null) {
							synchronized (statsIndex) {
								if (statsIndex[vi] == null) {
									statsIndex[vi] = new CountsPerTaxid(node.getLevel(), node.getTaxId(), initialReadSize);
								}
								stats = statsIndex[vi];
							}
						}
						synchronized (stats) {
							stats.reads++;
							stats.readsKmers += readKmers;
							stats.readsBPs += entry.readSize;
							stats.errorSum += err;
							stats.errorSquaredSum += err * err;
							stats.classErrorSum += classErr;
							stats.classErrorSquaredSum += classErr * classErr;
						}
					} else if (getLogger().isWarnEnabled()) {
						getLogger().warn("Missing database entry for tax node: " + node);
					}
				}
			}
		}

		return found;
	}
}
