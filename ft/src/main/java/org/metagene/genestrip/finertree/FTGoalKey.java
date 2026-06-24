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
    @MDDescription("Counts the *k*-mers per taxid directly from the underlying genomic files given a corresponding *k*-mer is in the database at all.")
    DB_QUALITY_COUNTS("dbqualcounts", false),
    @MDDescription("Generate a CSV file with intrinsic quality metrics per tax id for a project's database. "
            + "For each tax id, *k*-mers are counted directly from the underlying genomic files and compared with the *k*-mers stored in the database. "
            + "The resulting CSV file contains precision and recall values: "
            + "*precision* measures how specifically the stored *k*-mers are assigned to the correct tax id (true positives / (true positives + false positives)); "
            + "*recall* measures how completely the genomic *k*-mers of a tax id are represented in the database (true positives / (true positives + false negatives)). "
            + "The result is stored under `<base dir>/projects/<project_name>/csv`.")
    DB_QUALITY("dbquality", false);

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
