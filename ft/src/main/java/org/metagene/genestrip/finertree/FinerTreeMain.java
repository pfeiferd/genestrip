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

public abstract class FinerTreeMain<P extends FTProject> extends Main<P> {
    @Override
    protected FinerTreeMaker<P> createMaker(P project) {
        return new FinerTreeMaker<P>(project);
    }

    public static void main(String[] args) {
        new FinerTreeMain<FTProject>() {
            @Override
            protected FTProject createProject(GSCommon config, String name, String key, String[] fastqFiles, String csvFile, File csvDir, File fastqResDir, String taxids, Properties commandLineProps, GSGoalKey forGoal, String dbPath, boolean quietInit) {
                return new FTProject(config, name, key, fastqFiles, csvFile, csvDir, fastqResDir, taxids, commandLineProps, forGoal, dbPath, quietInit);
            }
        }.parseAndRun(args);
    }
}
