package org.metagene.genestrip.gen;

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
	private Project project;
	private Config config;
	private String target;
	private String[] restArgs;
	private Generator generator;

	public Main() {
		options = createOptions();
	}

	public void parse(String[] args) throws ParseException {
		try {
			CommandLine line = new DefaultParser().parse(options, args);

			String baseDir = line.getOptionValue("d", "./data");
			config = new Config(new File(baseDir));

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
			if (fastqName != null) {
				fastqFile = new File(fastqName);
			}

			restArgs = line.getArgs();
			project = new Project(config, restArgs[0], q, k, db, fastqFile);
			generator = new Generator(getProject());
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public String getTarget() {
		return target;
	}

	public Project getProject() {
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

		Option fastq = Option.builder("f").hasArg().argName("fastq file")
				.desc("Fastq file in case of filtering or classfication (regarding goals 'filter' and 'classify').")
				.build();
		options.addOption(fastq);

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
				case "cleanall":
					generator.cleanAll(restArgs[i]);
					break;
				case "clean":
					generator.clean(restArgs[i]);
					break;
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

	public Generator getGenerator() {
		return generator;
	}

	public static void main(String[] args) {
		new Main().parseAndRun(args);
	}
}
