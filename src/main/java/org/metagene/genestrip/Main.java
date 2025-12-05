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
import java.util.Properties;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.metagene.genestrip.make.GoalKey;

public class Main {
	private final Options options;
	private GSProject project;
	private String target;
	private String[] restArgs;
	private GSMaker maker;

	public Main() {
		options = createOptions();
	}

	public void parse(String[] args) throws ParseException, IOException {
		CommandLine line = new DefaultParser().parse(options, args);

		if (line.hasOption("v")) {
			Package p = Main.class.getPackage();
			System.out.println(p.getImplementationTitle() + " version " + p.getImplementationVersion());
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

		restArgs = line.getArgs();
		if (restArgs.length == 0) {
			throw new ParseException("Missing project name.");
		}
		String projectName = restArgs[0];

		Properties props = line.getOptionProperties("C");

		project = new GSProject(new GSCommon(new File(baseDir)), projectName, key, fastqFiles, mapFilePath, resFolder,
				resFolder, taxids, props, null, dbPath, false);
		project.setDownloadFastqs(download);
		project.setDownloadFastqsToCommon(downloadToCommon);

		maker = createMaker(project);
	}

	protected GSMaker createMaker(GSProject project) {
		return new GSMaker(project);
	}

	public String getTarget() {
		return target;
	}

	public GSProject getProject() {
		return project;
	}

	public Options getOptions() {
		return options;
	}

	public String[] getRestArgs() {
		return restArgs;
	}

	public String getDefaultGoal() {
		return "show";
	}

	protected Options createOptions() {
		Options options = new Options();

		Option version = Option.builder("v").hasArg(false).desc("Print version.").build();
		options.addOption(version);

		Option baseDir = Option.builder("d").hasArg().argName("base dir")
				.desc("Base directory for all data files. The default is './data'.").build();
		options.addOption(baseDir);

		Option target = Option.builder("t").hasArg().argName("target")
				.desc("Generation target ('make', 'clean' or 'cleanall'). The default is 'make'.").build();
		options.addOption(target);

		Option fastqs = Option.builder("f").hasArgs().valueSeparator(',').argName("fqfile1,fqfile2,...").desc(
				"Input fastq files as paths or URLs to be processed via the goals 'filter', 'match' or 'matchlr'. When a URL is given, the fastq file will not be downloaded but data streaming will be applied unless '-l' or '-ll' is given.")
				.build();
		options.addOption(fastqs);

		Option key = Option.builder("k").hasArgs().argName("key")
				.desc("Key used as a prefix for naming result files in conjuntion with '-f'.").build();
		options.addOption(key);

		Option csv = Option.builder("m").hasArg().argName("fqmap").desc(
				"Mapping file with a list of fastq files to be processed via the goals 'filter', 'match' or 'matchlr'. Each line of the file must have the format '<key> <URL or path to fastq file>'.")
				.build();
		options.addOption(csv);

		Option resFolder = Option.builder("r").hasArg().argName("path").desc(
				"Common store folder for filtered fastq files and result files created via the goals 'filter', 'match' or 'matchlr'. The defaults are '<base dir>/projects/<project name>/fastq' and '<base dir>/projects/<project name>/csv', respectively.")
				.build();
		options.addOption(resFolder);

		Option tax = Option.builder("tx").hasArg().argName("taxids").desc(
				"List of tax ids separated by ',' (but no blanks) for the goal 'db2fastq'. A tax id may have the suffix '+', which means that taxonomic descendants from the project's database will be included.")
				.build();
		options.addOption(tax);

		Option download = Option.builder("l").hasArg(false).argName("load").desc(
				"Download fastqs from URLs to '<base dir>/projects/<project name>/fastq' instead of streaming them for the goals 'filter', 'match' and 'matchlr'.")
				.build();
		options.addOption(download);

		Option downloadToCommon = Option.builder("ll").hasArg(false).argName("loadToCommon").desc(
				"Download fastqs from URLs to '<base dir>/fastq' instead of streaming them for the goals 'filter', 'match' and 'matchlr'.")
				.build();
		options.addOption(downloadToCommon);

		Option db = Option.builder("db").hasArg().argName("database").desc(
				"Path to filtering or matching database for the goals 'filter' or 'match', 'matchlr', 'dbinfo' and 'db2fastq' for use without project context.")
				.build();
		options.addOption(db);

		Option propsOption = Option.builder("C").hasArgs().valueSeparator('=').argName("key>=<value")
				.desc("To set Genestrip configuration paramaters via the command line.").build();
		options.addOption(propsOption);

		return options;
	}

	public void parseAndRun(String[] args) {
		try {
			parse(args);
			run(System.err);
		} catch (IOException | ParseException e) {
			HelpFormatter formatter = new HelpFormatter();
			Package p = Main.class.getPackage();			
			formatter.printHelp(p.getImplementationTitle() + " [options] <project> [<goal1> <goal2>...]", getOptions());
			System.out.flush();
			System.out.println();
			e.printStackTrace();
		} finally {
			if (maker != null) {
				maker.dumpAll();
			}
		}
	}

	public void run(PrintStream out) {
		String[] restArgs = getRestArgs();
		if (restArgs.length == 1) {
			restArgs = new String[] { null, maker.getDefaultGoalKey().getName() };
		}
		for (int i = 1; i < restArgs.length; i++) {
			GoalKey goalKey = maker.getKeyByName(restArgs[i]);
			if (goalKey != null) {
				out.println("Executing target " + getTarget() + " for goal " + goalKey.getName() + "...");
				switch (getTarget()) {
				case "cleantotal":
					maker.cleanTotal(goalKey);
					break;
				case "cleanall":
					maker.cleanAll(goalKey);
					break;
				case "clean":
					maker.clean(goalKey);
					break;
				default:
				case "make":
					maker.make(goalKey);
					break;
				}
				out.println("Done with target " + getTarget() + " for goal " + restArgs[i]);
			} else {
				out.println("Omitting unknown goal " + restArgs[i]);
			}
		}
	}

	public GSMaker getMaker() {
		return maker;
	}

	public static void main(String[] args) {
		new Main().parseAndRun(args);
	}
}
