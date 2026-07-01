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

import org.metagene.genestrip.*;

import java.io.File;
import java.util.Properties;

/**
 * Command line entry point for the finer-tree (FT) extension. Behaves like {@link Main} but builds a
 * {@link FinerTreeMaker} so the FT-specific goals are available.
 *
 * @param <P> the concrete FT project type
 */
public abstract class FinerTreeMain<P extends FTProject> extends Main<P> {
    /**
     * Creates the finer-tree command-line launcher.
     */
    protected FinerTreeMain() {
    }

    @Override
    protected FinerTreeMaker<P> createMaker(P project) {
        return new FinerTreeMaker<P>(project);
    }

    /**
     * Runs the finer-tree extension on the default {@link FTProject} type from the command line.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        new FinerTreeMain<FTProject>() {
            @Override
            protected FTProject createProject(GSCommon config, String name, String key, String[] fastqFiles, String csvFile, File csvDir, File fastqResDir, String taxids, Properties commandLineProps, GSGoalKey forGoal, String dbPath, boolean quietInit) {
                return new FTProject(config, name, key, fastqFiles, csvFile, csvDir, fastqResDir, taxids, commandLineProps, forGoal, dbPath, quietInit);
            }
        }.parseAndRun(args);
    }
}
