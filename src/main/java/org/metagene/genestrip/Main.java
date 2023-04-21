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
import org.metagene.genestrip.tax.AssemblySummaryReader.FTPEntryQuality;

public class Main {
	private final Options options;
	private GSProject project;
	private GSConfig config;
	private String target;
	private String[] restArgs;
	private GSMaker generator;

	public Main() {
		options = createOptions();
	}

	public void parse(String[] args) throws ParseException {
		try {
			CommandLine line = new DefaultParser().parse(options, args);

			String baseDir = line.getOptionValue("d", "./data");
			config = new GSConfig(new File(baseDir));

			String db = line.getOptionValue("db", null);
			String kmer = line.getOptionValue("k", "0");
			int k = Integer.valueOf(kmer);

			String qs = line.getOptionValue("q", null);
			FTPEntryQuality q = null;
			if (qs != null) {
				q = FTPEntryQuality.valueOf(qs.toUpperCase());
			}

			target = line.getOptionValue("t", "make");

			String fastqName = line.getOptionValue("f");
			File fastqFile = null;
			if (fastqName != null && !fastqName.trim().isEmpty()) {
				fastqFile = new File(fastqName.trim());
			}

			boolean ignoreKrakenOutFilter = line.hasOption('i');

			restArgs = line.getArgs();
			String projectName = restArgs[0];
			File resFolder = null;

			String resStr = line.getOptionValue("r");
			if (resStr != null) {
				resFolder = new File(resStr);
			}

			project = new GSProject(config, projectName, q, k, db, fastqFile, resFolder, !ignoreKrakenOutFilter);
			generator = new GSMaker(getProject());
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
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
				.desc("Base directrory of files, default is './data'.").build();
		options.addOption(baseDir);

		Option dbPath = Option.builder("db").hasArg().argName("path")
				.desc("Path to Kraken database folder, default is from 'Config.properties'.").build();
		options.addOption(dbPath);

		Option kMerSize = Option.builder("k").hasArg().argName("size")
				.desc("K-mer size, default is from 'Config.properties'.").build();
		options.addOption(kMerSize);

		Option quality = Option.builder("q").hasArg().argName("quality").desc(
				"Quality of fasta files from NCBI ('COMLPETE_LATEST', 'COMPLETE', 'LATEST', 'NONE'), default is from 'Config.properties'.")
				.build();
		options.addOption(quality);

		Option target = Option.builder("t").hasArg().argName("target")
				.desc("Generation target ('make', 'clean', 'cleanall'), default is 'make'.").build();
		options.addOption(target);

		Option fastq = Option.builder("f").hasArg().argName("fqfile")
				.desc("Fastq file in case of filtering or classfication (regarding goals 'filter' and 'classify').")
				.build();
		options.addOption(fastq);

		Option ignoreKrakenOutFilter = Option.builder("i").argName("ignoretaxids").desc(
				"Whether to ignore taxids from project when running goal 'krakenrescount'. If set only taxids from project will be counted in result file, default is 'don't ignore'.")
				.build();
		options.addOption(ignoreKrakenOutFilter);

		Option resFolder = Option.builder("r").hasArg().argName("path").desc(
				"Store folder files created via goals 'filter' and 'classify', default is '<base directory >/projects/<project name>/res'.")
				.build();
		options.addOption(resFolder);

		return options;
	}

	public void parseAndRun(String[] args) {
		try {
			parse(args);
			run(System.err);
		} catch (ParseException e) {
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
			restArgs = new String[] { null, generator.getDefaultGoalName() };
		}
		Set<String> goalNames = generator.getGoalNames();
		for (int i = 1; i < restArgs.length; i++) {
			if (goalNames.contains(restArgs[i])) {
				out.println("Executing target " + getTarget() + " for goal " + restArgs[i]);
				switch (getTarget()) {
				case "cleanTotal":
					generator.cleanTotal(restArgs[i]);
					break;
				case "cleanall":
					generator.cleanAll(restArgs[i]);
					break;
				case "clean":
					generator.clean(restArgs[i]);
					break;
				default:
				case "make":
					generator.make(restArgs[i]);
					break;
				}
				out.println("Done with target " + getTarget() + " for goal " + restArgs[i]);
			} else {
				out.println("Omitting unknown goal " + restArgs[i]);
			}
		}
	}

	public GSMaker getGenerator() {
		return generator;
	}

	public static void main(String[] args) {
		new Main().parseAndRun(args);
	}
}
