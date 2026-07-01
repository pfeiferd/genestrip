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
package org.metagene.genestrip;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Properties;

import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.make.ConfigKey;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.Project;
import org.metagene.genestrip.util.GSLogFactory;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;

/**
 * Represents a single Genestrip project: its configuration (loaded from command line properties and
 * {@code config.properties} files), its directory layout under the base directory, and its fastq/fasta
 * inputs. Also resolves output file names and locations for the various goals.
 */
public class GSProject extends Project {
    /** Database property key for the RefSeq release the database was built from. */
    public static final String REFSEQ_RELEASE = "refseq.release";
    /** Database property key for the Genestrip version that created the database. */
    public static final String GENESTRIP_VERSION = "genestrip.creationVersion";
    /** Database property key for the Genestrip title that created the database. */
    public static final String GENESTRIP_TITLE = "genestrip.creationTitle";
    /** Database property key for the database creation date. */
    public static final String DB_CREATION_DATE = "dbCreationDate";
    /** Database property key for the MD5 checksum of the database. */
    public static final String DB_MD5 = "dbMD5";

    /** Name of the configuration properties file loaded from project and base directories. */
    public static final String CONFIG_PROPERTIES = "config.properties";

    /**
     * A category of project file that carries a file name suffix.
     */
    public interface FileType {
        /**
         * Returns the file name suffix associated with this file type.
         *
         * @return the file name suffix associated with this file type
         */
        public String getSuffix();
    }

    /**
     * The built-in Genestrip file types and their file name suffixes.
     */
    public enum GSFileType implements FileType {
        /** A filtered/result fastq file. */
        FASTQ_RES(".fastq"),
        /** A fastq input file. */
        FASTQ(".fastq"),
        /** A fasta input file. */
        FASTA(".fasta"),
        /** A CSV result file. */
        CSV(".csv"),
        /** A Kraken output file. */
        KRAKEN_OUT(".out"),
        /** A Kraken result output file. */
        KRAKEN_OUT_RES(".out"),
        /** A serialized object file. */
        SER(".ser"),
        /** A database archive file. */
        DB(".zip"),
        /** A serialized filter file. */
        FILTER(".ser"),
        /** A log file. */
        LOG(".log"),
        /** An SVG graphics file. */
        SVG(".svg");

        private final String suffix;

        private GSFileType(String suffix) {
            this.suffix = suffix;
        }

        @Override
        public String getSuffix() {
            return suffix;
        }
    }

    /**
     * Returns the Genestrip runtime version read from the jar manifest.
     *
     * @return the implementation version of the Genestrip jar, or {@code null} if not packaged
     */
    public static String getGenestripRuntimeVersion() {
        return GSProject.class.getPackage().getImplementationVersion();
    }

    /**
     * Returns the Genestrip runtime title read from the jar manifest.
     *
     * @return the implementation title of the Genestrip jar, or {@code null} if not packaged
     */
    public static String getGenestripRuntimeTitle() {
        return GSProject.class.getPackage().getImplementationTitle();
    }

    private final GSCommon common;
    private final String key;
    private final String[] fastqResources;
    private final String fastqMapFile;
    private final File csvDir;
    private final File fastqResDir;
    private final Properties[] properties;
    private final String dbPath;

    private boolean downloadFastqs;
    private boolean downloadFastqsToCommon;

    /**
     * Creates a project with the given system configuration and project name, using defaults for all
     * other settings.
     *
     * @param config the shared system configuration
     * @param name   the project name (also the project directory name)
     */
    public GSProject(GSCommon config, String name) {
        this(config, name, null, null, null, null, null, null, null, null, null, false);
    }

    /**
     * Creates a project with the given system configuration and project name, optionally suppressing
     * informational logging during initialization.
     *
     * @param config the shared system configuration
     * @param name   the project name (also the project directory name)
     * @param quiet  if {@code true}, suppresses informational logging during initialization
     */
    public GSProject(GSCommon config, String name, boolean quiet) {
        this(config, name, null, null, null, null, null, null, null, null, null, quiet);
    }

    /**
     * Creates a project with the given result key and fastq inputs.
     *
     * @param config     the shared system configuration
     * @param name       the project name (also the project directory name)
     * @param key        key used as a prefix for result file names
     * @param fastqFiles fastq/fasta input paths or URLs, may be {@code null}
     */
    public GSProject(GSCommon config, String name, String key, String[] fastqFiles) {
        this(config, name, key, fastqFiles, null, null, null, null, null, null, null, false);
    }

    /**
     * Creates a project with the given result key and fastq inputs, optionally suppressing
     * informational logging during initialization.
     *
     * @param config     the shared system configuration
     * @param name       the project name (also the project directory name)
     * @param key        key used as a prefix for result file names
     * @param fastqFiles fastq/fasta input paths or URLs, may be {@code null}
     * @param quiet      if {@code true}, suppresses informational logging during initialization
     */
    public GSProject(GSCommon config, String name, String key, String[] fastqFiles, boolean quiet) {
        this(config, name, key, fastqFiles, null, null, null, null, null, null, null, quiet);
    }

    /**
     * Creates a project whose fastq/fasta inputs are listed in the given mapping (CSV) file.
     *
     * @param config  the shared system configuration
     * @param name    the project name (also the project directory name)
     * @param csvFile mapping file listing fastq/fasta inputs, may be {@code null}
     */
    public GSProject(GSCommon config, String name, String csvFile) {
        this(config, name, null, null, csvFile, null, null, null, null, null, null, false);
    }

    /**
     * Creates a fully configured project. Command line properties as well as the
     * {@code config.properties} files in the project and base directories are loaded and validated,
     * and the logger is configured accordingly.
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
    public GSProject(GSCommon config, String name, String key, String[] fastqFiles, String csvFile, File csvDir,
                     File fastqResDir, String taxids, Properties commandLineProps, GSGoalKey forGoal,
                     String dbPath, boolean quietInit) {
        super(name);
        this.common = config;
        this.fastqResDir = fastqResDir;
        this.key = key;
        this.fastqResources = fastqFiles;
        this.fastqMapFile = csvFile;
        this.csvDir = csvDir != null ? csvDir : new File(getProjectDir(), "csv");
        this.dbPath = dbPath;

        if (commandLineProps != null) {
            initConfigParams(commandLineProps); // Gotta do this here already, e.g. to impact logging early on.
            configureLogger();
        }

        properties = new Properties[3];
        properties[0] = commandLineProps == null ? new Properties() : commandLineProps;
        if (taxids != null) {
            properties[0].setProperty(GSConfigKey.TAX_IDS.getName(), taxids);
        }
        properties[1] = loadConfigProperties(getProjectDir(), quietInit);
        properties[2] = loadConfigProperties(getCommon().getBaseDir(), quietInit);

        initConfigParams(properties);
        configureLogger();

        for (Properties props : properties) {
            checkConfigProperties(props, forGoal);
        }
        if (!quietInit) {
            logParamMap();
        }
        GSProgressBarCreator.setGlobalUpdateInterval(intConfigValue(GSConfigKey.PROGRESS_BAR_UPDATE));
    }

    /**
     * Loads {@code config.properties} from the given directory, returning empty properties if the
     * file is missing or unreadable.
     *
     * @param dir   the directory to load {@code config.properties} from
     * @param quiet if {@code true}, suppresses info/warn logging about the file
     * @return the loaded properties, or empty properties if the file could not be read
     */
    protected Properties loadConfigProperties(File dir, boolean quiet) {
        Properties properties = new Properties();
        File configFile = new File(dir, CONFIG_PROPERTIES);
        if (!quiet && getLogger().isInfoEnabled()) {
            getLogger().info("Loading config file '" + configFile + "'.");
        }
        try (InputStream is = new FileInputStream(configFile)) {
            properties.load(is);
        } catch (IOException e) {
            if (!quiet && getLogger().isWarnEnabled()) {
                getLogger().warn("Could not read config file '" + configFile + "'.");
            }
        }
        return properties;
    }

    /**
     * Configures the logger's log level from the {@code LOG_LEVEL} configuration value.
     */
    protected void configureLogger() {
        GSLogFactory.getInstance().setLogLevel(stringConfigValue(GSConfigKey.LOG_LEVEL));
    }

    @Override
    protected ConfigKey[] getConfigKeys() {
        return GSConfigKey.values();
    }

    /**
     * Returns the key used as a prefix for result file names.
     *
     * @return the key used as a prefix for result file names, may be {@code null}
     */
    public String getKey() {
        return key;
    }

    /**
     * Returns the explicit database path for project-less operation.
     *
     * @return the explicit database path for project-less operation, may be {@code null}
     */
    public String getDBPath() {
        return dbPath;
    }

    /**
     * Returns the configured fastq/fasta input resources.
     *
     * @return the configured fastq/fasta input resources, may be {@code null}
     */
    public String[] getFastqResources() {
        return fastqResources;
    }

    /**
     * Hook for subclasses to supply additional input resources; returns {@code null} by default.
     *
     * @return the extra input resources, or {@code null} if none
     */
    // For override...
    public StreamingResourceStream getExtraResources() {
        return null;
    }

    /**
     * Hook for subclasses to supply the key for the extra resources; returns {@code null} by default.
     *
     * @return the key for the extra resources, or {@code null} if none
     */
    // For override...
    public String getExtraResourcesKey() {
        return null;
    }

    /**
     * Returns the mapping file listing fastq/fasta inputs.
     *
     * @return the mapping file listing fastq/fasta inputs, may be {@code null}
     */
    public String getFastqMapFile() {
        return fastqMapFile;
    }

    /**
     * Returns whether fastq inputs should be downloaded.
     *
     * @return {@code true} if fastq inputs should be downloaded
     */
    public boolean isDownloadFastqs() {
        return downloadFastqs;
    }

    /**
     * Returns whether downloaded fastq inputs should be stored in the common directory.
     *
     * @return {@code true} if downloaded fastq inputs should be stored in the common directory
     */
    public boolean isDownloadFastqsToCommon() {
        return downloadFastqsToCommon;
    }

    /**
     * Sets whether fastq inputs should be downloaded.
     *
     * @param downloadFastqs {@code true} to download fastq inputs
     */
    public void setDownloadFastqs(boolean downloadFastqs) {
        this.downloadFastqs = downloadFastqs;
    }

    /**
     * Sets whether downloaded fastq inputs should be stored in the common directory.
     *
     * @param downloadFastqsToCommon {@code true} to store downloaded fastq inputs in the common directory
     */
    public void setDownloadFastqsToCommon(boolean downloadFastqsToCommon) {
        this.downloadFastqsToCommon = downloadFastqsToCommon;
    }

    /**
     * Resolves the storage directory for files of the given type.
     *
     * @param type the file type to resolve the directory for
     * @return the directory in which files of the given type are stored
     * @throws IllegalArgumentException if the type is not a recognized {@link GSFileType}
     */
    public File getDirForType(FileType type) {
        if (type instanceof GSFileType) {
            GSFileType gsFileType = (GSFileType) type;
            switch (gsFileType) {
                case FASTQ_RES:
                    return getFastqResDir();
                case FASTQ:
                    return getFastqDir();
                case FASTA:
                    return getFastaDir();
                case CSV:
                case SVG:
                    return getResultsDir();
                case KRAKEN_OUT:
                case KRAKEN_OUT_RES:
                    return getKrakenOutDir();
                case SER:
                case DB:
                case FILTER:
                    return getDBDir();
                case LOG:
                    return getLogDir();
            }
        }
        throw new IllegalArgumentException("Illegal FileType: " + type);
    }

    /**
     * Resolves the gzip-compressed output file for the given goal and file type.
     *
     * @param goal the goal the file belongs to
     * @param type the file type
     * @return the resolved output file
     */
    public File getOutputFile(String goal, FileType type) {
        return getOutputFile(goal, type, true);
    }

    /**
     * Resolves the output file for the given goal and file type.
     *
     * @param goal the goal the file belongs to
     * @param type the file type
     * @param gzip if {@code true}, the file name gets a {@code .gz} extension
     * @return the resolved output file
     */
    public File getOutputFile(String goal, FileType type, boolean gzip) {
        return getOutputFile(goal, null, null, type, gzip);
    }

    /**
     * Resolves the gzip-compressed output file for the given goal, base file name and file type.
     *
     * @param goal     the goal the file belongs to
     * @param baseFile the base file name, may be {@code null}
     * @param type     the file type
     * @return the resolved output file
     */
    public File getOutputFile(String goal, String baseFile, FileType type) {
        return getOutputFile(goal, null, baseFile, type, true);
    }

    /**
     * Resolves the output file for the given goal, key, base file name and file type, in the
     * directory associated with the file type.
     *
     * @param goal     the goal the file belongs to
     * @param key      the result key, may be {@code null}
     * @param baseFile the base file name, may be {@code null}
     * @param type     the file type
     * @param gzip     if {@code true}, the file name gets a {@code .gz} extension
     * @return the resolved output file
     */
    public File getOutputFile(String goal, String key, String baseFile, FileType type, boolean gzip) {
        return getOutputFile(getDirForType(type), goal, key, baseFile, type, gzip);
    }

    /**
     * Resolves the output file for the given goal, key, base file name and file type in the given
     * directory, using the default project prefix.
     *
     * @param dir      the directory to place the file in
     * @param goal     the goal the file belongs to
     * @param key      the result key, may be {@code null}
     * @param baseFile the base file name, may be {@code null}
     * @param type     the file type
     * @param gzip     if {@code true}, the file name gets a {@code .gz} extension
     * @return the resolved output file
     */
    public File getOutputFile(File dir, String goal, String key, String baseFile, FileType type, boolean gzip) {
        return getOutputFile(dir, getOutputFilePrefix(goal), goal, key, baseFile, type, gzip);
    }

    /**
     * Builds an output file in the given directory whose name is composed from an optional project
     * prefix, a goal/key infix and the base file name (with the input suffix stripped), followed by
     * the type's suffix and an optional {@code .gz} extension.
     *
     * @param dir      the directory to place the file in
     * @param project  the project prefix, may be {@code null}
     * @param goal     the goal the file belongs to
     * @param key      the result key, may be {@code null}
     * @param baseFile the base file name, may be {@code null}
     * @param type     the file type
     * @param gzip     if {@code true}, the file name gets a {@code .gz} extension
     * @return the resolved output file
     */
    public File getOutputFile(File dir, String project, String goal, String key, String baseFile, FileType type,
                              boolean gzip) {
        String baseName = baseFile == null ? "" : getFileBaseName(baseFile);
        if (baseName.startsWith(getName() + "_")) {
            baseName = baseName.substring(getName().length() + 1);
        }
        String infix = getOutputFileGoalPrefix(goal, key);
        if (!infix.isEmpty()) {
            if (!baseName.isEmpty()) {
                baseName = infix + "_" + baseName;
            } else {
                baseName = infix;
            }
        }
        return new File(dir, (project == null ? "" : project) + baseName + type.getSuffix() + (gzip ? ".gz" : ""));
    }

    /**
     * Strips any compression and known file type suffix from the given file name.
     *
     * @param projectFileName the file name to strip
     * @return the given file name stripped of any {@code .gz}/{@code .gzip} extension and any known
     *         {@link GSFileType} suffix
     */
    public String getFileBaseName(String projectFileName) {
        String baseName = projectFileName;
        if (baseName.endsWith(".gz")) {
            baseName = baseName.substring(0, baseName.length() - 3);
        } else if (baseName.endsWith(".gzip")) {
            baseName = baseName.substring(0, baseName.length() - 5);
        }
        for (GSFileType ft : GSFileType.values()) {
            String suffix = ft.getSuffix();
            if (baseName.endsWith(suffix)) {
                baseName = baseName.substring(0, baseName.length() - suffix.length());
            }
        }
        return baseName;
    }

    /**
     * Builds the project-specific prefix for output file names.
     *
     * @param goal the goal the file belongs to
     * @return the output file name prefix
     */
    protected String getOutputFilePrefix(String goal) {
        return getName() + "_";
    }

    /**
     * Builds the goal/key infix used in output file names, URL-encoding the key (truncated to 256
     * characters) and joining it to the goal with an underscore when both are present.
     *
     * @param goal the goal name, may be {@code null}
     * @param key  the result key, may be {@code null}
     * @return the goal/key infix for output file names
     */
    protected String getOutputFileGoalPrefix(String goal, String key) {
        if (key != null) {
            try {
                key = URLEncoder.encode(key, "UTF-8");
                if (key.length() > 256) {
                    key = key.substring(0, 256);
                }
            } catch (UnsupportedEncodingException e) {
                throw new RuntimeException(e);
            }
        }
        if (goal == null) {
            return key == null ? "" : key;
        } else {
            return key == null ? goal : goal + "_" + key;
        }
    }

    /**
     * Returns the shared system configuration.
     *
     * @return the shared system configuration
     */
    public GSCommon getCommon() {
        return common;
    }

    /**
     * Returns the base {@code projects} directory containing all project directories.
     *
     * @return the base {@code projects} directory containing all project directories
     */
    public File getProjectsDir() {
        return new File(getCommon().getBaseDir(), "projects");
    }

    /**
     * Returns this project's directory.
     *
     * @return this project's directory
     */
    public File getProjectDir() {
        return new File(getProjectsDir(), getName());
    }

    /**
     * Returns the directory holding this project's fasta files.
     *
     * @return the directory holding this project's fasta files
     */
    public File getFastaDir() {
        return new File(getProjectDir(), "fasta");
    }

    /**
     * Returns the directory holding this project's fastq files.
     *
     * @return the directory holding this project's fastq files
     */
    public File getFastqDir() {
        return new File(getProjectDir(), "fastq");
    }

    /**
     * Returns the directory for filtered/result fastq files.
     *
     * @return the directory for filtered/result fastq files, defaulting to the fastq directory
     */
    public File getFastqResDir() {
        return fastqResDir == null ? getFastqDir() : fastqResDir;
    }

    /**
     * Returns the directory holding this project's database files.
     *
     * @return the directory holding this project's database files
     */
    public File getDBDir() {
        return new File(getProjectDir(), "db");
    }

    /**
     * Returns the database archive file.
     *
     * @return the database archive file
     */
    public File getDBFile() {
        return getOutputFile(GSGoalKey.DB.getName(), GSFileType.DB, false);
    }

    /**
     * Returns the database info CSV file.
     *
     * @return the database info CSV file
     */
    public File getDBInfoFile() {
        return getOutputFile(GSGoalKey.DBINFO.getName(), GSFileType.CSV, false);
    }

    /**
     * Returns the temporary database info CSV file.
     *
     * @return the temporary database info CSV file
     */
    public File getTempDBInfoFile() {
        return getOutputFile(GSGoalKey.TEMP_DBINFO.getName(), GSFileType.CSV, false);
    }

    /**
     * Returns the directory holding Kraken output files.
     *
     * @return the directory holding Kraken output files
     */
    public File getKrakenOutDir() {
        return new File(getProjectDir(), "krakenout");
    }

    /**
     * Returns the {@code taxids.txt} file listing the project's tax ids.
     *
     * @return the {@code taxids.txt} file listing the project's tax ids
     */
    public File getTaxIdsFile() {
        return new File(getProjectDir(), "taxids.txt");
    }

    /**
     * Returns the {@code additional.txt} file listing additional tax ids.
     *
     * @return the {@code additional.txt} file listing additional tax ids
     */
    public File getAdditionalFile() {
        return new File(getProjectDir(), "additional.txt");
    }

    /**
     * Returns the {@code categories.txt} file listing tax id categories.
     *
     * @return the {@code categories.txt} file listing tax id categories
     */
    public File getCategoriesFile() {
        return new File(getProjectDir(), "categories.txt");
    }

    /**
     * Returns the directory holding CSV/result files.
     *
     * @return the directory holding CSV/result files
     */
    public File getResultsDir() {
        return csvDir;
    }

    /**
     * Returns the directory holding log files.
     *
     * @return the directory holding log files
     */
    public File getLogDir() {
        return new File(getProjectDir(), "log");
    }

    /**
     * Resolves the filter file for the given goal.
     *
     * @param <P>  the project type
     * @param goal the goal whose key names the filter file
     * @return the resolved filter file
     */
    public <P extends GSProject> File getFilterFile(Goal<P> goal) {
        return getOutputFile(goal.getKey().getName(), GSFileType.FILTER);
    }

    /**
     * Resolves a fasta file, trying the path as-is, then relative to the project's fasta directory,
     * then relative to the common fasta directory.
     *
     * @param fastaFilePath the fasta file path to resolve
     * @return the first existing file, or {@code null} if none is found
     */
    public File fastaFileFromPath(String fastaFilePath) {
        if (fastaFilePath != null) {
            File fasta = new File(fastaFilePath);
            if (fasta.exists()) {
                return fasta;
            }
            fasta = new File(getFastaDir(), fastaFilePath);
            if (fasta.exists()) {
                return fasta;
            }
            fasta = new File(getCommon().getFastaDir(), fastaFilePath);
            if (fasta.exists()) {
                return fasta;
            }
        }
        return null;
    }

    /**
     * Resolves a fastq file, trying the path as-is, then relative to the project's fastq directory,
     * then relative to the common fastq directory.
     *
     * @param fastqFilePath the fastq file path to resolve
     * @return the first existing file, or {@code null} if none is found
     */
    public File fastqFileFromPath(String fastqFilePath) {
        if (fastqFilePath != null) {
            File fastq = new File(fastqFilePath);
            if (fastq.exists()) {
                return fastq;
            }
            fastq = new File(getFastqDir(), fastqFilePath);
            if (fastq.exists()) {
                return fastq;
            }
            fastq = new File(getCommon().getFastqDir(), fastqFilePath);
            if (fastq.exists()) {
                return fastq;
            }
        }
        return null;
    }

    /**
     * Resolves fastq (or fasta) inputs from a path that may be a direct file or a glob pattern,
     * searching the given path, the project directory and the common directory in turn.
     *
     * @param fastqFilePath the file path or glob pattern to resolve
     * @param fastqType if {@code true}, resolves against fastq directories, otherwise fasta directories
     * @return the matching files, or {@code null} if none are found
     */
    public List<File> fastqFilesFromPath(String fastqFilePath, boolean fastqType) {
        if (fastqFilePath != null) {
            File fastq = new File(fastqFilePath);
            if (fastq.exists()) {
                return Collections.singletonList(fastq);
            }
            fastq = new File(fastqType ? getFastqDir() : getFastaDir(), fastqFilePath);
            if (fastq.exists()) {
                return Collections.singletonList(fastq);
            }
            fastq = new File(fastqType ? getCommon().getFastqDir() : getCommon().getFastaDir(), fastqFilePath);
            if (fastq.exists()) {
                return Collections.singletonList(fastq);
            }
            List<File> list = null;
            list = findFilesByGlobPattern(null, fastqFilePath);
            if (list != null) {
                return list;
            }
            list = findFilesByGlobPattern(fastqType ? getFastqDir() : getFastaDir(), fastqFilePath);
            if (list != null) {
                return list;
            }
            list = findFilesByGlobPattern(fastqType ? getCommon().getFastqDir() : getCommon().getFastaDir(), fastqFilePath);
            if (list != null) {
                return list;
            }
        }
        return null;
    }

    /**
     * Lists the files in the pattern's parent directory (optionally under {@code rootDir}) whose
     * names match the glob part of the pattern.
     *
     * @param rootDir the root directory to resolve the pattern against, or {@code null} to use the pattern's own path
     * @param pattern the glob pattern whose parent directory is listed and whose name part is matched
     * @return the matching files, or {@code null} if none match
     */
    protected List<File> findFilesByGlobPattern(File rootDir, String pattern) {
        List<File> res = null;
        File dir;
        if (rootDir == null) {
            dir = new File(pattern).getParentFile();
        } else {
            dir = new File(rootDir, pattern).getParentFile();
        }
        if (dir != null) {
            File[] files = dir.listFiles();
            if (files != null && files.length > 0) {
                pattern = new File(pattern).getName();
                PathMatcher matcher = FileSystems.getDefault().getPathMatcher("glob:" + pattern);
                for (File file : files) {
                    Path path = file.toPath();
                    if (matcher.matches(path.getFileName())) {
                        if (res == null) {
                            res = new ArrayList<File>();
                        }
                        res.add(file);
                    }
                }
            }
        }
        return res;
    }
}
