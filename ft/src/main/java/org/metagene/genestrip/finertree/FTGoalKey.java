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

public enum FTGoalKey implements GoalKey {
    @MDDescription("Create additional folders in `<base dir>/<project>` like `tex`.")
    FTSETUP("ftsetup"),
    @MDDescription("Generate the dendrograms from *k*-mer intersection counts using agglomerative clustering.")
    DENDROGRAM("dendrogram"),
    @MDDescription("Generate LaTeX extracts for depicting the dendrograms from `dendrogram`.")
    DENDRO_LATEX("dendrolatex", true),
    @MDDescription("Store which *k*-mers belongs to which species for all *k*-mers under the rank genus (and potentially other ranks depending on configuration) in a Bloom filter.")
    KMER_INDEX_BLOOM("kmerindexbloom"),
    @MDDescription("Count the number of joint *k*-mers between any two species per genus rank (and potentially other ranks depending on configuration).")
    INTERSECT_COUNT("intersectcount"),
    @MDDescription("Save the number of joint *k*-mers between any two species per genus rank (and potentially other ranks depending on configuration) to CSV files along with resulting Jaccard-indices.")
    INTERSECT_CSV("intersectcsv", true),
    @MDDescription("Load the Bloom filter computed via the goal `kmerindexbloom`.")
    LOAD_KMER_INDEX("loadkmerindex"),
    @MDDescription("Store the Bloom filter computed via the goal `kmerindexbloom`.")
    STORE_KMER_INDEX("storekmerindex"),
    @MDDescription("Update the database by integrating the refined taxonomy tree and reassigning *k*-mers under the genus ranks accordingly.")
    UPDATE_STORE_GOAL("updatestore"),
    @MDDescription("Store the updated database.")
    FTDB("ftdb", true),
    @MDDescription("Write information on the updated database content to a CSV file.")
    FTDBINFO("ftdbinfo", true),
    @MDDescription("Load the updated database.")
    LOAD_FTDB("loadftdb"),
    @MDDescription("Merge a project's LaTeX extracts from `dendrolatex` into one LaTeX document.")
    ALLINONE_LATEX("allinonelatex", true),
    @MDDescription("Analyze fastq files according to Genestrip's `matchres` but with a Genestrip-FT database instead.")
    FTMATCHRES("ftmatchres"),
    @MDDescription("Analyze fastq files according to Genestrip's `match` but with a Genestrip-FT database instead.")
    FTMATCH("ftmatch", true),
    @MDDescription("Generate fastq files according to Genestrip's `db2fastq` but from a Genestrip-FT database instead.")
    FTDB2FASTQ("ftdb2fastq", true),
    @MDDescription("Same as goal `clear`, but also clears `tex` the folder.")
    FTCLEAR("ftclear", true),
    @MDDescription("Same as `svgtaxtree` but for an FT database.")
    FT_SVG_TAX_TREE("ftsvgtaxtree", true),
    @MDDescription("Counts the *k*-mers per taxid directly from the underlying genomic files given a corresponding *k*-mer is in the database at all.")
    KMERS_PER_TAXID("tax2kmers", false),
    @MDDescription("Same as `tax2kmers` but for FT database.")
    FT_KMERS_PER_TAXID("fttax2kmers", false),
    @MDDescription("TODO")
    FT_CONSISTENCY("ftcons", false),
    @MDDescription("TODO")
    DB_CONSISTENCY("dbcons", false);

    private final boolean forUser;
    private final String name;

    private FTGoalKey(String name) {
        this(name, false);
    }

    private FTGoalKey(String name, boolean forUser) {
        this.name = name;
        this.forUser = forUser;
    }

    @Override
    public boolean isTransClean() {
        return true;
    }

    public boolean isForUser() {
        return forUser;
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public String toString() {
        return name;
    }

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
