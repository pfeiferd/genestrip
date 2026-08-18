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

import org.junit.Assume;
import org.junit.Test;
import org.metagene.genestrip.*;
import org.metagene.genestrip.make.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Collection;

/**
 * Documentation-generation helper implemented as JUnit tests so it can be run on demand (typically
 * only for releases). It regenerates the Markdown tables of FT configuration parameters and goals as
 * well as the Graphviz goal-dependency graph of the finer-tree goal set.
 */
// Rendered as a test because normally not needed - just for releases...
public class GenDocFiles {
	/**
	 * Regenerates the {@code ConfigParams.md} and {@code Goals.md} documentation files from
	 * {@link FTConfigKey} and {@link FTGoalKey} in the project directory.
	 */
	@Test
	public void writeDokuFiles() {
		File configParamsFile = new File(getProjectDir(), "ConfigParams.md");
		try (PrintStream ps = new PrintStream(configParamsFile)) {
			FTConfigKey.printMDConfigParamInfo(ps, null);
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}

		File goalInfoFile = new File(getProjectDir(), "Goals.md");
		try (PrintStream ps = new PrintStream(goalInfoFile)) {
			FTGoalKey.printGoalInfo(ps);
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Regenerates the {@code GoalGraph.gv.txt} Graphviz file describing the finer-tree goal graph,
	 * emitting one node per relevant goal and directed edges for their dependencies (dotted for weak
	 * dependencies), and renders it as {@code GoalGraph.svg} via {@link #writeSVGFile(File)}.
	 *
	 * @throws IOException if writing the graph file or rendering it fails
	 */
	@Test
	public void writeGraphFile() throws IOException {
		File graphFile = new File(getProjectDir(), "GoalGraph.gv.txt");
		try (PrintStream ps = new PrintStream(graphFile)) {
			GSCommon config = new GSCommon(getBaseDir());
			FTProject project = new FTProject(config, "virus", null, null, null, null, null, null, null, null, null, false);
			FinerTreeMaker maker = new FinerTreeMaker(project);

			ps.println("digraph regexp {");
			ps.println("rankdir=\"BT\"");
			ps.println("fontname=\"Helvetica,Arial,sans-serif\"");
			ps.println("node [fontname=\"Helvetica,Arial,sans-serif\", style=rounded, shape=box]");
			//ps.println("edge [fontname=\"Helvetica,Arial,sans-serif\"]");
			Collection<Goal<GSProject>> goals = maker.getGoals();
			for (Goal<GSProject> goal : goals) {
				if (!isRelevantGoal(goal, goals)) {
					continue;
				}
				if (goal.getKey().equals(GSGoalKey.SETUP)) {
					continue;
				}
				ps.print(goal.getKey().getName());
				ps.print(" [label=\"");
				if (goal instanceof ObjectGoal) {
					ps.print("o:");
				}
				else if (goal instanceof FileDownloadGoal) {
					ps.print("d:");					
				}
				else if (goal instanceof FileGoal) {
					ps.print("f:");					
				}
				ps.print(goal.getKey().getName());
				ps.print("\"");
				GoalKey key = goal.getKey();
				if (key instanceof GSGoalKey) {
					ps.print(" style=\"rounded, dashed\"");
				}
				else if (key instanceof FTGoalKey) {
					if (((FTGoalKey) key).isForUser()) {
						ps.print(" style=\"bold\"");
					}
				}
				ps.println("];");
			}

			for (Goal<GSProject> goal : goals) {
				if (!isRelevantGoal(goal, goals)) {
					continue;
				}
				if (goal.getKey() instanceof GSGoalKey) {
					continue;
				}
				if (goal.getKey().equals(GSGoalKey.SETUP)) {
					continue;
				}
				for (Goal<GSProject> to : goal.getDependencies()) {
					if (to.getKey().equals(GSGoalKey.SETUP)) {
						continue;
					}
					ps.print(goal.getKey().getName());
					ps.print(" -> ");
					ps.print(to.getKey().getName());
					if (goal.isWeakDependency(to)) {
						ps.print(" [style=dotted]");						
					}
					ps.println(";");
				}
			}
			ps.println("}");
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}
		writeSVGFile(graphFile);
	}

	/**
	 * Renders the given Graphviz file next to itself as an SVG file (with the {@code .gv...} part of its
	 * name replaced by {@code .svg}) by running {@code dot} from Graphviz (https://graphviz.org/).
	 * If {@code dot} is not on the PATH, the SVG file is left as it is and the calling test is skipped
	 * rather than failed.
	 *
	 * @param graphFile the Graphviz file to render
	 * @throws IOException if {@code dot} could not be run to completion or reported a failure
	 */
	public static void writeSVGFile(File graphFile) throws IOException {
		String name = graphFile.getName();
		int index = name.indexOf(".gv");
		File svgFile = new File(graphFile.getParentFile(), (index < 0 ? name : name.substring(0, index)) + ".svg");

		ProcessBuilder builder = new ProcessBuilder("dot", "-Tsvg", graphFile.getAbsolutePath(), "-o",
				svgFile.getAbsolutePath());
		builder.redirectErrorStream(true);
		Process process;
		try {
			process = builder.start();
		} catch (IOException e) {
			String message = "Graphviz's 'dot' is not on the PATH, so " + svgFile
					+ " was not regenerated. Install Graphviz (https://graphviz.org/) for that.";
			System.err.println(message);
			Assume.assumeNoException(message, e);
			return;
		}
		StringBuilder output = new StringBuilder();
		try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				output.append(line).append('\n');
			}
		}
		int exitValue;
		try {
			exitValue = process.waitFor();
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new IOException("Interrupted while waiting for 'dot' to create " + svgFile, e);
		}
		if (exitValue != 0) {
			throw new IOException("'dot' failed with exit value " + exitValue + " when creating " + svgFile + ": "
					+ output);
		}
	}

	private boolean isRelevantGoal(Goal goal, Collection<Goal<GSProject>> goals) {
		if (goal.getKey() instanceof FTGoalKey) {
			return true;
		}
		else if (goal.getKey() instanceof GSGoalKey) {
			for (Goal<GSProject> g : goals) {
				if (g.getKey() instanceof FTGoalKey) {
					for (Goal<GSProject> g2 : g.getDependencies()) {
						if (g2 == goal) {
							return true;
						}
					}
				}
			}
			return false;
		}
		else {
			return true;
		}
	}

	/**
	 * Returns the base directory holding the release's project data, resolved as the {@code data}
	 * subdirectory of the project directory.
	 *
	 * @return the base data directory
	 */
	// To find the base directory for any project file provided as part of the
	// release.
	// (This may have to be done differently for you own projects.)
	public static File getBaseDir() {
		return new File(getProjectDir(), "data");
	}

	/**
	 * Returns the project directory, taken from the {@code project.directory} system property if set
	 * and otherwise derived from the location of this class's compiled code.
	 *
	 * @return the project directory
	 */
	public static File getProjectDir() {
		String projectDir = System.getProperty("project.directory");
		if (projectDir != null) {
			return new File(projectDir);
		}
		String relPath = GenDocFiles.class.getProtectionDomain().getCodeSource().getLocation().getFile();
		return new File(relPath).getParentFile().getParentFile();
	}
}
