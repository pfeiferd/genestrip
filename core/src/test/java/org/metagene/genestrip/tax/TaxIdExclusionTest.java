package org.metagene.genestrip.tax;

import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;

import static org.junit.Assert.assertEquals;

/**
 * Pins what a {@code taxids.txt} means, in particular the exclusion a line prefixed with {@code -}
 * expresses.
 * <p>
 * The exclusion has been in {@link TaxIdCollector} from early on but was documented nowhere, so it
 * was neither used nor guarded; the {@code streptonamed} database of the FT paper is the first thing
 * that depends on it. What it has to keep doing is stated here rather than in prose: an exclusion
 * removes a whole branch however deep it reaches, it works regardless of where in the file it
 * stands, and it does nothing at all when it names a branch no inclusion brought in.
 * <p>
 * The tree is written as a pair of miniature NCBI dumps rather than mocked, so that the real parser
 * of {@link TaxTree} is exercised alongside the collector: a taxonomy read differently would break
 * the exclusion just as surely as a collector that forgot it.
 */
public class TaxIdExclusionTest {
    @Rule
    public final TemporaryFolder folder = new TemporaryFolder();

    private TaxTree tree;
    private TaxIdCollector collector;

    /**
     * A miniature Streptococcus: the genus with one named species carrying a strain, beside the two
     * unnamed buckets that the real genus has and that {@code streptonamed} leaves out.
     *
     * <pre>
     * 1 root
     *  +- 1301 genus Streptococcus
     *      +- 1313 species S. pneumoniae
     *      |    +- 170187 strain
     *      +- 2608887 no rank "unclassified Streptococcus"
     *      |    +- 1306 species S. sp.
     *      +- 83426 no rank "environmental samples"
     *      |    +- 9999 species S. sp. env.
     *      +- 1314 species S. pyogenes
     * </pre>
     */
    @Before
    public void setUp() throws IOException {
        File dir = folder.newFolder("taxonomy");
        write(new File(dir, "nodes.dmp"),
                "1\t|\t1\t|\tno rank\t|",
                "1301\t|\t1\t|\tgenus\t|",
                "1313\t|\t1301\t|\tspecies\t|",
                "170187\t|\t1313\t|\tstrain\t|",
                "2608887\t|\t1301\t|\tno rank\t|",
                "1306\t|\t2608887\t|\tspecies\t|",
                "83426\t|\t1301\t|\tno rank\t|",
                "9999\t|\t83426\t|\tspecies\t|",
                "1314\t|\t1301\t|\tspecies\t|");
        write(new File(dir, "names.dmp"),
                "1\t|\troot\t|\t\t|\tscientific name\t|",
                "1301\t|\tStreptococcus\t|\t\t|\tscientific name\t|",
                "1313\t|\tStreptococcus pneumoniae\t|\t\t|\tscientific name\t|",
                "2608887\t|\tunclassified Streptococcus\t|\t\t|\tscientific name\t|",
                "83426\t|\tenvironmental samples\t|\t\t|\tscientific name\t|");
        tree = new TaxTree(dir, false);
        collector = new TaxIdCollector(tree);
    }

    private static void write(File file, String... lines) throws IOException {
        try (PrintWriter pw = new PrintWriter(file, StandardCharsets.UTF_8.name())) {
            for (String line : lines) {
                pw.println(line);
            }
        }
    }

    /** Runs a {@code taxids.txt} through the same two steps {@code TaxNodesGoal} runs it through. */
    private Set<String> resolve(String... lines) throws IOException {
        File file = folder.newFile();
        write(file, lines);
        Set<TaxIdNode> excludes = new HashSet<>();
        Set<TaxIdNode> included = collector.readFromFile(file, excludes);
        Set<String> res = new TreeSet<>();
        for (TaxIdNode node : collector.completeAndExclude(included, excludes, null)) {
            res.add(node.getTaxId());
        }
        return res;
    }

    private static Set<String> taxIds(String... ids) {
        return new TreeSet<>(java.util.Arrays.asList(ids));
    }

    /** Without an exclusion the genus brings everything below it, buckets included. */
    @Test
    public void testAGenusBringsItsWholeSubtree() throws IOException {
        assertEquals(taxIds("1301", "1313", "170187", "2608887", "1306", "83426", "9999", "1314"),
                resolve("1301"));
    }

    /** The case streptonamed is built on: the genus without its two unnamed branches. */
    @Test
    public void testAnExclusionRemovesABranchOfAnIncludedTaxon() throws IOException {
        assertEquals(taxIds("1301", "1313", "170187", "1314"),
                resolve("1301", "-2608887", "-83426"));
    }

    /** The branch goes in full: excluding the bucket takes the species under it as well. */
    @Test
    public void testAnExclusionReachesAllTheWayDown() throws IOException {
        assertEquals(taxIds("1301", "1313", "170187", "83426", "9999", "1314"),
                resolve("1301", "-2608887"));
    }

    /** Inclusions are expanded first and exclusions subtracted after, so the order cannot matter. */
    @Test
    public void testOrderOfLinesDoesNotMatter() throws IOException {
        assertEquals(resolve("1301", "-2608887"), resolve("-2608887", "1301"));
    }

    /** An exclusion may name a taxon nothing brought in; it then simply has no effect. */
    @Test
    public void testAnUncoveredExclusionDoesNothing() throws IOException {
        assertEquals(taxIds("1313", "170187"), resolve("1313", "-2608887"));
    }

    /** Excluding what was included leaves nothing, rather than the ancestor of it. */
    @Test
    public void testExcludingTheIncludedTaxonItselfLeavesNothing() throws IOException {
        assertEquals(taxIds(), resolve("1313", "-1313"));
    }

    /** Comments are ignored, whole-line and trailing alike, which is how our own files are written. */
    @Test
    public void testCommentsAreIgnored() throws IOException {
        assertEquals(taxIds("1301", "1313", "170187", "1314"),
                resolve("# the genus, without what has no name",
                        "1301        # everything below Streptococcus ...",
                        "-2608887    # ... but not the unnamed drafts",
                        "-83426      # ... nor the environmental samples",
                        ""));
    }

    /** A tab-separated two-column list is read as it stands, the tax id being the last field. */
    @Test
    public void testTaxIdIsReadAfterTheLastTab() throws IOException {
        assertEquals(taxIds("1313", "170187"), resolve("Streptococcus pneumoniae\t1313"));
    }
}
