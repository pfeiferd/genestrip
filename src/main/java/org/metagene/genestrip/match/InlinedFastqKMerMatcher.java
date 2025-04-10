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
package org.metagene.genestrip.match;

import org.metagene.genestrip.ExecutionContext;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;
import org.metagene.genestrip.util.CGAT;

public class InlinedFastqKMerMatcher extends FastqKMerMatcher {
	public InlinedFastqKMerMatcher(KMerSortedArray<SmallTaxIdNode> kmerStore, int initialReadSize, int maxQueueSize, ExecutionContext bundle, boolean withProbs, int maxKmerResCounts, SmallTaxTree taxTree, int maxPaths, double maxReadTaxErrorCount, double maxReadClassErrorCount) {
		super(kmerStore, initialReadSize, maxQueueSize, bundle, withProbs, maxKmerResCounts, taxTree, maxPaths, maxReadTaxErrorCount, maxReadClassErrorCount);
	}
/*
	// Manually inlined stuff via Intellij - seems to help a bit.
	protected boolean matchRead(final MyReadEntry entry, final int index, final boolean reverse);
	{
		boolean found = false;
		int prints = 0;
		int readTaxErrorCount = taxTree == null ? -1 : 0;

		SmallTaxIdNode taxIdNode;
		int max = entry.readSize - k + 1;
		SmallTaxIdNode lastTaxid = null;
		int contigLen = 0;
		CountsPerTaxid stats = null;

		//ByteArrayUtil.println(entry.read, System.out);
		long kmer = -1;
		for (int i = 0; i < max; i++) {
			if (kmer == -1) {
				int c;
				kmer = 0;
				if (reverse) {
					for (int i1 = i + k - 1; i1 >= i; i1--) {
						// Inlined: res = Long.rotateLeft(res, 2);
						kmer = (kmer << 2) | (kmer >>> -2);
						c = CGAT.CGAT_REVERSE_JUMP_TABLE[entry.read[i1]];
						if (c == -1) {
							kmer = -1;
							i = i1;
							break;
						}
						kmer += c;
					}
				} else {
					int max1 = i + k;
					for (int i1 = i; i1 < max1; i1++) {
						// Inlined: res = Long.rotateLeft(res, 2);
						kmer = (kmer << 2) | (kmer >>> -2);
						c = CGAT.CGAT_JUMP_TABLE[entry.read[i1]];
						if (c == -1) {
							kmer = -1;
							i = i1;
							break;
						}
						kmer += c;
					}
				}
			} else {
				if (reverse) {
					int c = CGAT.CGAT_REVERSE_JUMP_TABLE[entry.read[i + k - 1]];
					if (c == -1) {
						kmer = -1;
						i += k - 1;
					} else {
						kmer = (kmer >>> 2) | (((long) c) << CGAT.SHIFT_FILTERS_REVERSE[k]);
					}
				} else {
					int c = CGAT.CGAT_JUMP_TABLE[entry.read[i + k - 1]];
					if (c == -1) {
						kmer = -1;
						i += k - 1;
					} else {
						kmer = ((kmer << 2) & CGAT.SHIFT_FILTERS_STRAIGHT[k]) | (long) c;
					}
				}
			}
			if (kmer != -1) {
				taxIdNode = kmerStore.getLongInlined(kmer, entry.indexPos);
				if (readTaxErrorCount != -1) {
					if (taxIdNode == null) {
						readTaxErrorCount++;
						if (maxReadTaxErrorCount >= 0) {
							if ((maxReadTaxErrorCount >= 1 && readTaxErrorCount > maxReadTaxErrorCount)
									|| (readTaxErrorCount > maxReadTaxErrorCount * max)) {
								readTaxErrorCount = -1;
							}
						}
					} else {
						taxIdNode.incCount(index, entry.readNo, consumers);

						boolean found1 = false;
						for (int i1 = 0; i1 < entry.usedPaths; i1++) {
							if (taxTree.isAncestorOf(taxIdNode, entry.readTaxIdNode[i1])) {
								entry.readTaxIdNode[i1] = taxIdNode;
								found1 = true;
								break;
							} else if (taxTree.isAncestorOf(entry.readTaxIdNode[i1], taxIdNode)) {
								found1 = true;
								break;
							}
						}
						if (!found1) {
							if (entry.usedPaths < maxPaths) {
								entry.readTaxIdNode[entry.usedPaths] = taxIdNode;
								entry.usedPaths++;
							}
						}
					}
				}
				if (taxIdNode != lastTaxid) {
					if (contigLen > 0) {
						if (out != null) {
							printKrakenStyleOut(entry, lastTaxid, contigLen, prints++, reverse);
						}
						if (stats != null) {
							synchronized (stats) {
								stats.contigs++;
								if (contigLen > stats.maxContigLen) {
									stats.maxContigLen = contigLen;
									int j = 0;
									for (; j < entry.readDescriptor.length - 1 && entry.readDescriptor[j] != 0; j++) {
										stats.maxContigDescriptor[j] = entry.readDescriptor[j];
									}
									stats.maxContigDescriptor[j] = 0;
								}
							}
						}
						contigLen = 0;
					}
				}
				contigLen++;
				lastTaxid = taxIdNode;
				if (taxIdNode != null) {
					short vi = taxIdNode.storeIndex;
					stats = statsIndex[vi];
					if (stats == null) {
						synchronized (statsIndex) {
							if (statsIndex[vi] == null) {
								statsIndex[vi] = new CountsPerTaxid(taxIdNode.taxId, initialReadSize);
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
						uniqueCounter.putInlined(kmer, taxIdNode.taxId, entry.indexPos[0]);
					}
				} else {
					stats = null;
				}
			}
		}
		if (found) {
			if (contigLen > 0) {
				if (out != null) {
					printKrakenStyleOut(entry, lastTaxid, contigLen, prints, reverse);
				}
				if (stats != null) {
					synchronized (stats) {
						stats.contigs++;
						if (contigLen > stats.maxContigLen) {
							stats.maxContigLen = contigLen;
							for (int j = 0; j < entry.readDescriptorSize; j++) {
								stats.maxContigDescriptor[j] = entry.readDescriptor[j];
							}
							stats.maxContigDescriptor[entry.readDescriptorSize] = 0;
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
				SmallTaxIdNode node = entry.readTaxIdNode[0];
				for (int i = 1; i <= ties; i++) {
					node = taxTree.getLeastCommonAncestor(node, entry.readTaxIdNode[i]);
				}
				entry.classNode = node;
				// For 'readKmers', I decided to count in the k-mers from 'entry.readTaxIdNode[0]' and not just 'node'.
				// (They only differ in case of a tie anyways.) But if there is tie, then the k-mers from one of the tie's nodes
				// solidify the LCA in a sense - so the counts from one the involved paths are included.
				int readKmers = ties > 0 ? taxTree.sumCounts(entry.readTaxIdNode[0], index, entry.readNo) : entry.counts[0];
				int classErrC = max - readKmers;
				if (maxReadClassErrorCount < 0 || (maxReadClassErrorCount >= 1 && classErrC <= maxReadClassErrorCount)
						|| (classErrC <= maxReadClassErrorCount * max)) {
					double err = ((double) readTaxErrorCount) / max;
					double classErr = ((double) classErrC) / max;
					entry.classNode = node;
					short vi = node.getStoreIndex();
					stats = statsIndex[vi];
					if (stats == null) {
						synchronized (statsIndex) {
							if (statsIndex[vi] == null) {
								statsIndex[vi] = new CountsPerTaxid(node.getTaxId(), initialReadSize);
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
				}
			}
		}

		return found;
	}
 */
}
