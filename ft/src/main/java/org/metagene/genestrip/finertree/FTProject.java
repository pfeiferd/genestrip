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

/**
 * A {@link GSProject} for the finer-tree (FT) extension. Adds the FT-specific configuration keys and
 * the TeX file type and directory used when generating finer-tree output.
 */
public class FTProject extends GSProject {
    /**
     * FT-specific file types and their suffixes (currently the TeX output file).
     */
    public enum FTFileType implements FileType {
        /**
         * The TeX output file type.
         */
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

    /**
     * Creates a fully configured finer-tree project, delegating to the {@link GSProject} constructor.
     *
     * @param config           the shared system configuration
     * @param name             the project name (also the project directory name)
     * @param key              key used as a prefix for result file names
     * @param fastqFiles       fastq/fasta input paths or URLs, may be {@code null}
     * @param csvFile          mapping file listing fastq/fasta inputs, may be {@code null}
     * @param csvDir           output directory for CSV/result files; defaults to the project's csv folder
     * @param fastqResDir      output directory for filtered fastq files, may be {@code null}
     * @param taxids           comma-separated tax ids, applied as the {@code taxids} config value
     * @param commandLineProps configuration properties supplied on the command line, may be {@code null}
     * @param forGoal          the goal the project is being set up for, used when validating config keys
     * @param dbPath           explicit database path for project-less operation, may be {@code null}
     * @param quietInit        if {@code true}, suppresses informational logging during initialization
     */
    public FTProject(GSCommon config, String name, String key, String[] fastqFiles, String csvFile, File csvDir,
                     File fastqResDir, String taxids, Properties commandLineProps, GSGoalKey forGoal,
                     String dbPath, boolean quietInit) {
        super(config, name, key, fastqFiles, csvFile, csvDir, fastqResDir, taxids, commandLineProps, forGoal, dbPath, quietInit);
    }

    /**
     * @return the base Genestrip configuration keys plus the FT-specific {@link FTConfigKey} values
     */
    @Override
    protected ConfigKey[] getConfigKeys() {
        if (configKeys == null) {
            List<ConfigKey> keys = new ArrayList<ConfigKey>(Arrays.asList(super.getConfigKeys()));
            keys.addAll(Arrays.asList(FTConfigKey.values()));
            configKeys = keys.toArray(new ConfigKey[keys.size()]);
        }
        return configKeys;
    }

    /**
     * Returns the directory used for TeX output.
     *
     * @return the {@code tex} subdirectory of the project directory used for TeX output
     */
    public File getTeXDir() {
        return new File(getProjectDir(), "tex");
    }

    /**
     * Resolves the directory for FT-specific file types (TeX), delegating to the superclass otherwise.
     */
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
