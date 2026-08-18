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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Collection;

import org.junit.Assume;
import org.junit.Test;
import org.metagene.genestrip.make.FileDownloadGoal;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.ResultReporter;

// Rendered as a test because normally not needed - just for releases...
public class GenDocFiles {
	@Test
	public void writeDokuFiles() {
		File configParamsFile = new File(APITest.getProjectDir(), "ConfigParams.md");
		try (PrintStream ps = new PrintStream(configParamsFile)) {
			GSConfigKey.printMDConfigParamInfo(ps, null);
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}

		File goalInfoFile = new File(APITest.getProjectDir(), "Goals.md");
		try (PrintStream ps = new PrintStream(goalInfoFile)) {
			GSGoalKey.printGoalInfo(ps);
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}

		File csvColFile = new File(APITest.getProjectDir(), "CSVColumns.md");
		try (PrintStream ps = new PrintStream(csvColFile)) {
			ResultReporter.printMDColumnInfo(ps);
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}
	}

	@Test
	public void writeGraphFile() throws IOException {
		File graphFile = new File(APITest.getProjectDir(), "GoalGraph.gv.txt");
		try (PrintStream ps = new PrintStream(graphFile)) {
			GSCommon config = new GSCommon(APITest.getBaseDir());
			GSProject project = new GSProject(config, "human_virus", null, new String[0]);
			GSMaker maker = new GSMaker(project);

			ps.println("digraph regexp {");
			ps.println("rankdir=\"BT\"");
			ps.println("fontname=\"Helvetica,Arial,sans-serif\"");
			ps.println("node [fontname=\"Helvetica,Arial,sans-serif\", style=rounded, shape=box]");
			//ps.println("edge [fontname=\"Helvetica,Arial,sans-serif\"]");
			Collection<Goal<GSProject>> goals = maker.getGoals();
			for (Goal<GSProject> goal : goals) {
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
				if (((GSGoalKey) goal.getKey()).isForUser()) {
					ps.print(" style=\"bold\"");
				}
				ps.println("];");
			}

			for (Goal<GSProject> goal : goals) {
				if (goal.getKey().equals(GSGoalKey.SETUP)) {
					continue;
				}
				for (Goal<GSProject> to : goal.getDependencies()) {
					if (to == null || to.getKey().equals(GSGoalKey.SETUP)) {
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

	// Renders the Graphviz file next to itself as an SVG file via 'dot' from Graphviz
	// (https://graphviz.org/). If 'dot' is not installed, the SVG file is left as it is and the
	// test is skipped rather than failed.
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
}
