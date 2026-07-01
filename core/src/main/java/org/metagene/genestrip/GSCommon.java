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

/**
 * Holds system-wide configuration rooted at a base directory and derives the standard shared
 * directory locations (common, fastq, fasta, refseq and genbank) used across projects.
 */
public class GSCommon {
	private final File baseDir;

	/**
	 * Creates the system configuration rooted at the given base directory.
	 *
	 * @param baseDir the base directory under which shared directories are located
	 * @throws IOException if the base directory cannot be resolved
	 */
	public GSCommon(File baseDir) throws IOException {
		this.baseDir = baseDir;
	}

	/**
	 * Returns the base directory under which all shared directories are located.
	 *
	 * @return the base directory
	 */
	public File getBaseDir() {
		return baseDir;
	}


	/**
	 * Returns the shared {@code common} directory located under the base directory.
	 *
	 * @return the shared {@code common} directory
	 */
	public File getCommonDir() {
		return new File(baseDir, "common");
	}

	/**
	 * Returns the shared {@code fastq} directory located under the base directory.
	 *
	 * @return the shared {@code fastq} directory
	 */
	public File getFastqDir() {
		return new File(baseDir, "fastq");
	}

	/**
	 * Returns the shared {@code fasta} directory located under the {@code common} directory.
	 *
	 * @return the shared {@code fasta} directory
	 */
	public File getFastaDir() {
		return new File(getCommonDir(), "fasta");
	}

	/**
	 * Returns the shared {@code refseq} directory located under the {@code common} directory.
	 *
	 * @return the shared {@code refseq} directory
	 */
	public File getRefSeqDir() {
		return new File(getCommonDir(), "refseq");
	}

	/**
	 * Returns the shared {@code genbank} directory located under the {@code common} directory.
	 *
	 * @return the shared {@code genbank} directory
	 */
	public File getGenbankDir() {
		return new File(getCommonDir(), "genbank");
	}
}
