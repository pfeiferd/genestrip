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
package org.metagene.genestrip.refseq;

import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxTree;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;

/**
 * Tests {@link ReworkingStoreFastaReader#foldUp(TaxIdNode, Rank, Set)}, which decides at which rank
 * a genome is filed.
 * <p>
 * The rule this pins down is what makes {@code foldTaxaBelow=species} usable for the quality
 * measures: a genome whose accession resolves below the species is filed at the species instead, so
 * that the species -- and not whichever strain node the taxonomy happens to carry -- becomes the
 * data taxon a *k*-mer is counted against.
 */
public class FoldTaxaBelowTest {
    @Rule
    public final TemporaryFolder folder = new TemporaryFolder();

    private TaxTree tree;

    /**
     * A miniature Streptococcus: genus, one species, one strain below it, and a second species with
     * no strain node at all -- the asymmetry the folding exists to remove.
     */
    @Before
    public void setUp() throws IOException {
        File dir = folder.newFolder("taxonomy");
        write(new File(dir, "nodes.dmp"),
                "1\t|\t1\t|\tno rank\t|",
                "1301\t|\t1\t|\tgenus\t|",
                "1313\t|\t1301\t|\tspecies\t|",
                "170187\t|\t1313\t|\tstrain\t|",
                "1314\t|\t1301\t|\tspecies\t|");
        write(new File(dir, "names.dmp"),
                "1\t|\troot\t|\t\t|\tscientific name\t|",
                "1301\t|\tStreptococcus\t|\t\t|\tscientific name\t|",
                "1313\t|\tStreptococcus pneumoniae\t|\t\t|\tscientific name\t|",
                "1314\t|\tStreptococcus pyogenes\t|\t\t|\tscientific name\t|");
        tree = new TaxTree(dir, false);
    }

    /** A strain is filed at its species, which is the whole point of the fold. */
    @Test
    public void testStrainFoldsOntoItsSpecies() {
        assertSame(node("1313"), ReworkingStoreFastaReader.foldUp(node("170187"), Rank.SPECIES, null));
    }

    /** A genome already at the fold rank stays where it is. */
    @Test
    public void testSpeciesStaysAtItself() {
        assertSame(node("1313"), ReworkingStoreFastaReader.foldUp(node("1313"), Rank.SPECIES, null));
    }

    /** Above the fold rank nothing is folded: a genus is not filed at a species. */
    @Test
    public void testGenusIsLeftAlone() {
        assertSame(node("1301"), ReworkingStoreFastaReader.foldUp(node("1301"), Rank.SPECIES, null));
    }

    /** Without a rank to fold at, every node is returned unchanged. */
    @Test
    public void testNoFoldRankLeavesTheNode() {
        assertSame(node("170187"), ReworkingStoreFastaReader.foldUp(node("170187"), null, null));
    }

    /**
     * The species must be among the requested taxa. A database that asks for the strain but not for
     * the species would otherwise have its genome filed at a node it never requested.
     */
    @Test
    public void testFoldOnlyOntoRequestedTaxa() {
        Set<TaxIdNode> requested = new HashSet<>(Collections.singletonList(node("1301")));
        assertSame(node("170187"), ReworkingStoreFastaReader.foldUp(node("170187"), Rank.SPECIES, requested));
        requested.add(node("1313"));
        assertSame(node("1313"), ReworkingStoreFastaReader.foldUp(node("170187"), Rank.SPECIES, requested));
    }

    /** A species without a strain node below it is unaffected -- it is already filed correctly. */
    @Test
    public void testSpeciesWithoutStrainIsUnaffected() {
        assertEquals("1314", ReworkingStoreFastaReader.foldUp(node("1314"), Rank.SPECIES, null).getTaxId());
    }

    private TaxIdNode node(String taxId) {
        return tree.getNodeByTaxId(taxId);
    }

    private static void write(File file, String... lines) throws IOException {
        try (PrintWriter pw = new PrintWriter(file, StandardCharsets.UTF_8.name())) {
            for (String line : lines) {
                pw.println(line);
            }
        }
    }
}
