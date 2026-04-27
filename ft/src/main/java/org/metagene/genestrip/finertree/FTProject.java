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

import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.ConfigKey;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;

public class FTProject extends GSProject {
    public enum FTFileType implements FileType {
        TEX(".tex");

        private final String suffix;

        private FTFileType(String suffix) {
            this.suffix = suffix;
        }

        public String getSuffix() {
            return suffix;
        }
    }

    private ConfigKey[] configKeys;

    public FTProject(GSCommon config, String name, String key, String[] fastqFiles, String csvFile, File csvDir,
                     File fastqResDir, String taxids, Properties commandLineProps, GSGoalKey forGoal,
                     String dbPath, boolean quietInit) {
        super(config, name, key, fastqFiles, csvFile, csvDir, fastqResDir, taxids, commandLineProps, forGoal, dbPath, quietInit);
    }

    @Override
    protected ConfigKey[] getConfigKeys() {
        if (configKeys == null) {
            List<ConfigKey> keys = new ArrayList<ConfigKey>(Arrays.asList(super.getConfigKeys()));
            keys.addAll(Arrays.asList(FTConfigKey.values()));
            configKeys = keys.toArray(new ConfigKey[keys.size()]);
        }
        return configKeys;
    }

    public File getTeXDir() {
        return new File(getProjectDir(), "tex");
    }

    @Override
    public File getDirForType(FileType type) {
        if (type instanceof FTFileType) {
            FTFileType ftFileType = (FTFileType) type;
            switch (ftFileType) {
                case TEX:
                    return getTeXDir();
            }
        }
        return super.getDirForType(type);
    }
}
