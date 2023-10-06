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
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class Main {
	private final Options options;
	private GSProject project;
	private GSConfig config;
	private String target;
	private String[] restArgs;
	private GSMaker maker;

	public Main() {
		options = createOptions();
	}

	public void parse(String[] args) throws ParseException, IOException {
		CommandLine line = new DefaultParser().parse(options, args);

		String baseDir = line.getOptionValue("d", "./data");
		config = new GSConfig(new File(baseDir));

		target = line.getOptionValue("t", "make");

		String fastqName = line.getOptionValue("f");
		File fastqFile = null;
		if (fastqName != null && !fastqName.trim().isEmpty()) {
			fastqFile = new File(fastqName.trim());
		}

		File resFolder = null;
		String resStr = line.getOptionValue("r");
		if (resStr != null) {
			resFolder = new File(resStr);
		}

		restArgs = line.getArgs();
		if (restArgs.length == 0) {
			throw new ParseException("Missing project name.");
		}
		String projectName = restArgs[0];

		project = new GSProject(config, projectName, 0, null, fastqFile, resFolder, resFolder);
		maker = new GSMaker(project);
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
//		options.addOption(new Option("v", "verbose", true, "Verbose information."));

		Option baseDir = Option.builder("d").hasArg().argName("path")
				.desc("Base directory for all data files. The default is './data'.").build();
		options.addOption(baseDir);

		Option target = Option.builder("t").hasArg().argName("target")
				.desc("Generation target ('make', 'clean', 'cleanall'). The default is 'make'.").build();
		options.addOption(target);

		Option fastq = Option.builder("f").hasArg().argName("fqfile").desc(
				"Input fastq file in case of filtering or k-mer matching (regarding goals 'filter' and 'match') or "
						+ "csv file with a list of fastq files to be processed via 'multimatch'. Each line should have the format '<prefix> <path to fastq file>'.")
				.build();
		options.addOption(fastq);

		Option resFolder = Option.builder("r").hasArg().argName("path").desc(
				"Common store folder for filtered fastq files and csv files created via the goals 'filter' and 'match'. The default is '<base directory>/projects/<project name>/fastq' and '<base directory>/projects/<project name>/csv' respectively.")
				.build();
		options.addOption(resFolder);

		return options;
	}

	public void parseAndRun(String[] args) {
		try {
			parse(args);
			run(System.err);
		} catch (IOException | ParseException e) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("genestrip [options] <project> [<goal1> <goal2>...]", getOptions());
			System.out.flush();
			System.out.println();
			e.printStackTrace();
		}
	}

	public void run(PrintStream out) {
		String[] restArgs = getRestArgs();
		if (restArgs.length == 1) {
			restArgs = new String[] { null, maker.getDefaultGoalName() };
		}
		Set<String> goalNames = maker.getGoalNames();
		for (int i = 1; i < restArgs.length; i++) {
			if (goalNames.contains(restArgs[i])) {
				out.println("Executing target " + getTarget() + " for goal " + restArgs[i]);
				switch (getTarget()) {
				case "cleantotal":
					maker.cleanTotal(restArgs[i]);
					break;
				case "cleanall":
					maker.cleanAll(restArgs[i]);
					break;
				case "clean":
					maker.clean(restArgs[i]);
					break;
				default:
				case "make":
					maker.make(restArgs[i]);
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
