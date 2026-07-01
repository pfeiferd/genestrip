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
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.metagene.genestrip.make.GoalKey;

/**
 * Command line entry point for Genestrip. Parses the command line options, builds the project and its
 * {@link GSMaker}, and runs the requested target (make/clean variants) for the given goals.
 *
 * @param <P> the concrete project type created by this launcher
 */
public abstract class Main<P extends GSProject> {
    private final Options options;
    private P project;
    private String target;
    private String[] restArgs;
    private GSMaker<P> maker;
    private boolean isolateGoals;

    /**
     * Creates a new instance, building the command line option set.
     */
    public Main() {
        options = createOptions();
    }

    /**
     * Parses the command line arguments, creating the project (with its configuration) and the maker.
     *
     * @param args the command line arguments
     * @throws ParseException if the arguments are malformed or the project name is missing
     * @throws java.io.IOException if project creation requires I/O that fails
     */
    public void parse(String[] args) throws ParseException, IOException {
        CommandLine line = new DefaultParser().parse(options, args);

        if (line.hasOption("v")) {
            System.out.println(GSProject.getGenestripRuntimeTitle() + " version " + GSProject.getGenestripRuntimeVersion());
        }

        String baseDir = line.getOptionValue("d", "./data");

        target = line.getOptionValue("t", "make");

        String key = line.getOptionValue("k");
        String[] fastqFiles = line.getOptionValues("f");
        String mapFilePath = line.getOptionValue("m");

        File resFolder = null;
        String resStr = line.getOptionValue("r");
        if (resStr != null) {
            resFolder = new File(resStr);
        }

        String taxids = line.getOptionValue("tx");

        boolean download = line.hasOption("l");
        boolean downloadToCommon = line.hasOption("ll");

        String dbPath = line.getOptionValue("db");

        isolateGoals = line.hasOption("i");

        restArgs = line.getArgs();
        if (restArgs.length == 0) {
            throw new ParseException("Missing project name.");
        }
        String projectName = restArgs[0];

        Properties props = line.getOptionProperties("C");

        project = createProject(new GSCommon(new File(baseDir)), projectName, key, fastqFiles, mapFilePath, resFolder,
                resFolder, taxids, props, null, dbPath, false);
        project.setDownloadFastqs(download);
        project.setDownloadFastqsToCommon(downloadToCommon);

        maker = createMaker(project);
    }

    /**
     * Creates the concrete project instance; implemented by subclasses to supply their project type.
     *
     * @param config           the common Genestrip configuration
     * @param name             the project name
     * @param key              the key used as a prefix for result file names, or {@code null}
     * @param fastqFiles       the input fastq/fasta file paths or URLs, or {@code null}
     * @param csvFile          the mapping file path, or {@code null}
     * @param csvDir           the directory for CSV result files, or {@code null}
     * @param fastqResDir      the directory for filtered fastq/fasta files, or {@code null}
     * @param taxids           the comma-separated list of tax ids, or {@code null}
     * @param commandLineProps the configuration properties set on the command line
     * @param forGoal          the goal the project is created for, or {@code null}
     * @param dbPath           the path to a filtering or matching database, or {@code null}
     * @param quietInit        whether to suppress output during initialization
     * @return the newly created project
     */
    protected abstract P createProject(GSCommon config, String name, String key, String[] fastqFiles, String csvFile, File csvDir,
                                       File fastqResDir, String taxids, Properties commandLineProps, GSGoalKey forGoal,
                                       String dbPath, boolean quietInit);

    /**
     * Creates the maker used to build and run the goal graph; may be overridden to use a subclass.
     *
     * @param project the project to create the maker for
     * @return the newly created maker
     */
    protected GSMaker<P> createMaker(P project) {
        return new GSMaker<P>(project);
    }

    /**
     * Returns the parsed generation target.
     *
     * @return the target name
     */
    public String getTarget() {
        return target;
    }

    /**
     * Returns the project created during parsing.
     *
     * @return the project
     */
    public P getProject() {
        return project;
    }

    /**
     * Returns the command line option set.
     *
     * @return the options
     */
    public Options getOptions() {
        return options;
    }

    /**
     * Returns the non-option arguments (project name followed by goal names).
     *
     * @return the remaining arguments
     */
    public String[] getRestArgs() {
        return restArgs;
    }

    /**
     * Returns the name of the default goal used when none is given.
     *
     * @return the default goal name
     */
    public String getDefaultGoal() {
        return "show";
    }

    /**
     * Builds the Apache Commons CLI option set understood by Genestrip.
     *
     * @return the option set
     */
    protected Options createOptions() {
        Options options = new Options();

        Option version = Option.builder("v").hasArg(false).desc("Print version.").build();
        options.addOption(version);

        Option baseDir = Option.builder("d").hasArg().argName("base dir")
                .desc("Base directory for all data files. The default is './data'.").build();
        options.addOption(baseDir);

        Option target = Option.builder("t").hasArg().argName("target")
                .desc("Generation target ('make', 'clean', 'cleanall' or 'cleantotal'). The default is 'make'.").build();
        options.addOption(target);

        Option isolateGoals = Option.builder("i").hasArg(false).argName("isolate goals")
                .desc("If several goals are given, it runs them without reusing intermediate results of object goals between the given goals. Isolation allows for freeing memory resources held by intermediate results but may also cause recomputation. The default is no isolation.").build();
        options.addOption(isolateGoals);

        Option fastqs = Option.builder("f").hasArgs().valueSeparator(',').argName("fqfile1,fqfile2,...").desc(
                        "Input fastq/fasta files as paths or URLs to be processed via the goals 'filter', 'match' or 'matchlr'. When a URL is given, the fastq/fasta file will not be downloaded but data streaming will be applied unless '-l' or '-ll' is given.")
                .build();
        options.addOption(fastqs);

        Option key = Option.builder("k").hasArg().argName("key")
                .desc("Key used as a prefix for naming result files in conjuntion with '-f'.").build();
        options.addOption(key);

        Option csv = Option.builder("m").hasArg().argName("fqmap").desc(
                        "Mapping file with a list of fastq/fasta files to be processed via the goals 'filter', 'match' or 'matchlr'. Each line of the file must have the format '<key> <URL or path to fastq/fasta file>'.")
                .build();
        options.addOption(csv);

        Option resFolder = Option.builder("r").hasArg().argName("path").desc(
                        "Common store folder for filtered fastq/fasta files and result files created via the goals 'filter', 'match' or 'matchlr'. The defaults are '<base dir>/projects/<project name>/fastq' and '<base dir>/projects/<project name>/csv', respectively.")
                .build();
        options.addOption(resFolder);

        Option tax = Option.builder("tx").hasArg().argName("taxids").desc(
                        "List of tax ids separated by ',' (but no blanks) for the goal 'db2fastq/fasta'. A tax id may have the suffix '+', which means that taxonomic descendants from the project's database will be included. It can alternatively be set via the configuration parameter 'taxid'.")
                .build();
        options.addOption(tax);

        Option download = Option.builder("l").hasArg(false).argName("load").desc(
                        "Download fastq/fastas from URLs to '<base dir>/projects/<project name>/fastq' instead of streaming them for the goals 'filter', 'match' and 'matchlr'.")
                .build();
        options.addOption(download);

        Option downloadToCommon = Option.builder("ll").hasArg(false).argName("loadToCommon").desc(
                        "Download fastq/fastas from URLs to '<base dir>/fastq' instead of streaming them for the goals 'filter', 'match' and 'matchlr'.")
                .build();
        options.addOption(downloadToCommon);

        Option db = Option.builder("db").hasArg().argName("database").desc(
                        "Path to filtering or matching database for the goals 'filter' or 'match', 'matchlr', 'dbinfo' and 'db2fastq' for use without project context.")
                .build();
        options.addOption(db);

        Option propsOption = Option.builder("C").hasArgs().valueSeparator('=').argName("key>=<value")
                .desc("To set Genestrip configuration parameters via the command line.").build();
        options.addOption(propsOption);

        return options;
    }

    /**
     * Parses the arguments and runs the requested goals, printing the usage help on failure and
     * always disposing the maker's resources at the end.
     *
     * @param args the command line arguments
     */
    public void parseAndRun(String[] args) {
        try {
            parse(args);
            run(System.err);
        } catch (IOException | ParseException e) {
            HelpFormatter formatter = new HelpFormatter();
            Package p = Main.class.getPackage();
            formatter.printHelp("genestrip [options] <project> [<goal1> <goal2>...]", getOptions());
            System.out.flush();
            System.out.println();
            e.printStackTrace();
        } finally {
            if (maker != null) {
                maker.dumpAll();
            }
        }
    }

    /**
     * Resolves the goal names from the remaining arguments (falling back to the default goal) and
     * executes the parsed target ({@code make}, {@code clean}, {@code cleanall} or {@code cleantotal})
     * for them.
     *
     * @param out the stream for progress and status messages
     */
    public void run(PrintStream out) {
        String[] restArgs = getRestArgs();
        if (restArgs.length == 1) {
            restArgs = new String[]{null, maker.getDefaultGoalKey().getName()};
        }
        List<GoalKey> keys = new ArrayList<GoalKey>();
        for (int i = 1; i < restArgs.length; i++) {
            GoalKey goalKey = maker.getKeyByName(restArgs[i]);
            if (goalKey != null) {
                keys.add(goalKey);
            } else {
                out.println("Omitting unknown goal " + restArgs[i]);
            }
        }
        out.println("Executing target " + getTarget() + " for goals " + keys + "...");
        GoalKey[] ks = keys.toArray(new GoalKey[keys.size()]);
        switch (getTarget()) {
            case "cleantotal":
                maker.cleanTotal(ks);
                break;
            case "cleanall":
                maker.cleanAll(ks);
                break;
            case "clean":
                maker.clean(ks);
                break;
            default:
            case "make":
                maker.make(isolateGoals, ks);
                break;
        }
        out.println("Done with target " + getTarget() + " for goals " + keys);
    }

    /**
     * Returns the maker created during parsing.
     *
     * @return the maker, or {@code null} if parsing has not run
     */
    public GSMaker<P> getMaker() {
        return maker;
    }

    /**
     * Runs Genestrip on the default {@link GSProject} type from the command line.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        new Main<GSProject>() {
            protected GSProject createProject(GSCommon config, String name, String key, String[] fastqFiles, String csvFile, File csvDir,
                                              File fastqResDir, String taxids, Properties commandLineProps, GSGoalKey forGoal,
                                              String dbPath, boolean quietInit) {
                return new GSProject(config, name, key, fastqFiles, csvFile, csvDir,
                        fastqResDir, taxids, commandLineProps, forGoal,
                        dbPath, quietInit);
            }
        }.parseAndRun(args);
    }
}
