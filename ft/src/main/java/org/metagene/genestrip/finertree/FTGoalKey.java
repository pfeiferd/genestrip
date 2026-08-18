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
package org.metagene.genestrip.finertree;

import org.metagene.genestrip.make.MDDescription;
import org.metagene.genestrip.make.GoalKey;

import java.io.PrintStream;
import java.lang.annotation.Annotation;

/**
 * Enumeration of the goals added by the finer-tree (FT) extension. Constants annotated with
 * {@link MDDescription} contribute to the generated documentation.
 */
public enum FTGoalKey implements GoalKey {
    /** Generates the databases of Genestrip's {@code genall} and additionally the FT database. */
    @MDDescription("Generate the *k*-mer matching database and the filtering database according to Genestrip's `genall` and additionally the Genestrip-FT database with respect to the given project.")
    FTGENALL("ftgenall", true),
    /** Creates the additional project folders used by the FT extension, such as {@code tex}. */
    @MDDescription("Create additional folders in `<base dir>/<project>` like `tex`.")
    FTSETUP("ftsetup"),
    /** Generates the dendrograms from the k-mer intersection counts by agglomerative clustering. */
    @MDDescription("Generate the dendrograms from *k*-mer intersection counts using agglomerative clustering.")
    DENDROGRAM("dendrogram"),
    /** Generates the LaTeX extracts depicting the dendrograms of {@link #DENDROGRAM}. */
    @MDDescription("Generate LaTeX extracts for depicting the dendrograms from `dendrogram`.")
    DENDRO_LATEX("dendrolatex", true),
    /** Records in a Bloom filter which k-mer belongs to which species below the refinement ranks. */
    @MDDescription("Store which *k*-mers belongs to which species for all *k*-mers under the rank genus (and potentially other ranks depending on configuration) in a Bloom filter.")
    KMER_INDEX_BLOOM("kmerindexbloom"),
    /** Estimates how many entries the filter of {@link #KMER_INDEX_BLOOM} will have to hold. */
    @MDDescription("Estimate how many entries the `kmerindexbloom` filter will hold, by reading the sequences once and sketching the (*k*-mer, genome) pairs with HyperLogLog. Only made when `kmerIndexSizing` asks for that estimate.")
    KMER_INDEX_SIZE("kmerindexsize"),
    /** Counts the joint k-mers of any two species below a refinement rank. */
    @MDDescription("Count the number of joint *k*-mers between any two species per genus rank (and potentially other ranks depending on configuration).")
    INTERSECT_COUNT("intersectcount"),
    /** Writes the joint k-mer counts of {@link #INTERSECT_COUNT} and the Jaccard indices to CSV. */
    @MDDescription("Save the number of joint *k*-mers between any two species per genus rank (and potentially other ranks depending on configuration) to CSV files along with resulting Jaccard-indices.")
    INTERSECT_CSV("intersectcsv", true),
    /** Computes, per tree node with refined k-mers, a histogram over its k-mers' branching degrees. */
    @MDDescription("Compute, per tree node with refined *k*-mers, a histogram over the branching degrees of its *k*-mers, i.e. how many *k*-mers occur in exactly 1, 2, ... of the node's child subtrees (plus the trailing OTHER bucket as an additional child).")
    BRANCH_HISTO("branchhisto"),
    /** Writes the branching-degree histograms of {@link #BRANCH_HISTO} to CSV, one row per node. */
    @MDDescription("Write the branching-degree histograms from `branchhisto` to a CSV file, one row per tree node.")
    BRANCH_HISTO_CSV("branchhistocsv"),
    /** Aggregates the branching-degree histograms of {@link #BRANCH_HISTO} by rank into a CSV file. */
    @MDDescription("Aggregate the branching-degree histograms from `branchhisto` by taxonomic rank, one row per rank. Each node is turned into a relative branching-degree distribution (summing to one, excluding the OTHER column); for each branching degree from 1 up to 10, the mean, standard deviation, median and q1/q3 quartiles of that relative value across the rank's nodes are written as columns named `<degree>-avg`, `<degree>-stddev`, `<degree>-q1`, `<degree>-median`, `<degree>-q3`. Each row also starts with two summary blocks of the same five statistics: `childdeg-*` over the node's number of children in the database's taxonomy tree, and `kmerdeg-*` over the k-mer-weighted mean branching degree of the node (both excluding the OTHER column). Written to a CSV file.")
    BRANCH_HISTO_RANK_CSV("branchhistorankcsv"),
    /** Loads the Bloom filter computed by {@link #KMER_INDEX_BLOOM} from disk. */
    @MDDescription("Load the Bloom filter computed via the goal `kmerindexbloom`.")
    LOAD_KMER_INDEX("loadkmerindex"),
    /** Stores the Bloom filter computed by {@link #KMER_INDEX_BLOOM} to disk. */
    @MDDescription("Store the Bloom filter computed via the goal `kmerindexbloom`.")
    STORE_KMER_INDEX("storekmerindex"),
    /** Updates the database with the refined taxonomy tree and reassigns the k-mers accordingly. */
    @MDDescription("Update the database by integrating the refined taxonomy tree and reassigning *k*-mers under the genus ranks accordingly.")
    UPDATE_STORE_GOAL("ftupdatedb"),
    /** Stores the updated (FT) database. */
    @MDDescription("Store the updated database.")
    FTDB("ftdb", true),
    /** Writes information on the content of the updated (FT) database to a CSV file. */
    @MDDescription("Write information on the updated database content to a CSV file.")
    FTDBINFO("ftdbinfo", true),
    /** Loads the updated (FT) database. */
    @MDDescription("Load the updated database.")
    LOAD_FTDB("loadftdb"),
    /** Merges the LaTeX extracts of {@link #DENDRO_LATEX} into a single LaTeX document. */
    @MDDescription("Merge a project's LaTeX extracts from `dendrolatex` into one LaTeX document.")
    ALLINONE_LATEX("allinonelatex", true),
    /** Analyzes fastq files as Genestrip's {@code matchres} does, but against an FT database. */
    @MDDescription("Analyze fastq files according to Genestrip's `matchres` but with a Genestrip-FT database instead.")
    FTMATCHRES("ftmatchres"),
    /** Analyzes fastq files as Genestrip's {@code match} does, but against an FT database. */
    @MDDescription("Analyze fastq files according to Genestrip's `match` but with a Genestrip-FT database instead.")
    FTMATCH("ftmatch", true),
    /** Generates fastq files as Genestrip's {@code db2fastq} does, but from an FT database. */
    @MDDescription("Generate fastq files according to Genestrip's `db2fastq` but from a Genestrip-FT database instead.")
    FTDB2FASTQ("ftdb2fastq", true),
    /** Behaves like Genestrip's {@code clear} goal but also clears the {@code tex} folder. */
    @MDDescription("Same as goal `clear`, but also clears `tex` the folder.")
    FTCLEAR("ftclear", true),
    /** Behaves like Genestrip's {@code svgtaxtree} goal but for an FT database. */
    @MDDescription("Same as `svgtaxtree` but for an FT database.")
    FT_SVG_TAX_TREE("ftsvgtaxtree", true),
    /** Counts, per tax id, the genomic k-mers that are also in the ordinary database. */
    @MDDescription("Counts the *k*-mers per taxid directly from the underlying genomic files given a corresponding *k*-mer is in the database at all.")
    DB_QUALITY_COUNTS("dbqualcounts", false),
    /** Counts, per tax id, the genomic k-mers that are also in the FT database. */
    @MDDescription("Same as `dbqualcounts` but for an FT database.")
    FT_QUALITY_COUNTS("ftqualcounts", false),
    /** Writes the per-tax-id quality metrics of {@link #DB_QUALITY_COUNTS} to a CSV file. */
    @MDDescription("Write the per-taxid quality metrics (tp, tp+fp, tp+fn, precision and recall) derived from `dbqualcounts` to a CSV file.")
    DB_QUALITY("dbquality", false),
    /** Writes the per-tax-id quality metrics of {@link #FT_QUALITY_COUNTS} to a CSV file. */
    @MDDescription("Same as `dbquality` but for an FT database, i.e. derived from `ftqualcounts`.")
    FT_QUALITY("ftquality", false);

    private final boolean forUser;
    private final String name;

    private FTGoalKey(String name) {
        this(name, false);
    }

    private FTGoalKey(String name, boolean forUser) {
        this.name = name;
        this.forUser = forUser;
    }

    /**
     * Indicates whether this goal participates in transitive cleaning; FT goals always do.
     *
     * @return {@code true}, as FT goals are always transitively cleaned
     */
    @Override
    public boolean isTransClean() {
        return true;
    }

    /**
     * Indicates whether this goal is meant to be invoked directly by users.
     *
     * @return whether this goal is intended to be invoked directly by users
     */
    public boolean isForUser() {
        return forUser;
    }

    /**
     * Returns the textual name of this goal as used on the command line.
     *
     * @return the goal's name
     */
    @Override
    public String getName() {
        return name;
    }

    /**
     * Returns the goal's name.
     *
     * @return the goal's name
     */
    @Override
    public String toString() {
        return name;
    }

    /**
     * Prints a Markdown table listing the FT goals with their user flag and descriptions.
     *
     * @param ps the stream the Markdown table is written to
     */
    public static void printGoalInfo(PrintStream ps) {
        ps.print('|');
        ps.print("Name");
        ps.print('|');
        ps.print("User Goal");
        ps.print('|');
        ps.print("Description");
        ps.print('|');
        ps.println();

        ps.print('|');
        ps.print('-');
        ps.print('|');
        ps.print('-');
        ps.print('|');
        ps.print('-');
        ps.print('|');
        ps.println();

        for (FTGoalKey goalKey : FTGoalKey.values()) {
            ps.print('|');
            ps.print('`');
            ps.print(goalKey.getName());
            ps.print('`');
            ps.print('|');
            ps.print(goalKey.isForUser() ? "X" : "");
            ps.print('|');
            Annotation[] annotations;
            try {
                annotations = FTGoalKey.class.getField(goalKey.name()).getAnnotations();
            } catch (NoSuchFieldException e) {
                throw new RuntimeException(e);
            } catch (SecurityException e) {
                throw new RuntimeException(e);
            }
            for (Annotation annotation : annotations) {
                if (annotation instanceof MDDescription) {
                    ps.print(((MDDescription) annotation).value());
                    break;
                }
            }
            ps.print('|');
            ps.println();
        }
    }
}
