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
package org.metagene.genestrip.goals.refseq;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import java.io.File;

import org.junit.Test;
import org.metagene.genestrip.APITest;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSConfigKey.GenomesOnly;
import org.metagene.genestrip.GSProject;

/**
 * Tests the three-valued setting that says how finished an assembly must be to enter the database.
 */
public class GenomesOnlyConfigTest {

	private GSProject project() throws Exception {
		File baseDir = APITest.getBaseDir();
		return new GSProject(new GSCommon(baseDir), "human_virus", true);
	}

	/** Off by default: a database is not silently restricted to finished genomes. */
	@Test
	public void testDefaultIsOff() throws Exception {
		assertSame(GenomesOnly.OFF, project().configValue(GSConfigKey.GENOMES_ONLY));
	}

	/** Each of the three values must be settable and must arrive unchanged. */
	@Test
	public void testEachValueIsConfigurable() throws Exception {
		for (GenomesOnly value : GenomesOnly.values()) {
			GSProject project = project();
			assertTrue("must accept " + value, project.initConfigParam(GSConfigKey.GENOMES_ONLY, value));
			assertSame(value, project.configValue(GSConfigKey.GENOMES_ONLY));
		}
	}

	/**
	 * The names are what a user writes in a config.properties, so they are part of the interface and
	 * a rename has to be made here as well.
	 */
	@Test
	public void testConfigurationNames() {
		assertEquals("refseq.genomesOnly", GSConfigKey.GENOMES_ONLY.getName());
		assertEquals("off", GenomesOnly.OFF.getName());
		assertEquals("complete", GenomesOnly.COMPLETE.getName());
		assertEquals("chromosome", GenomesOnly.CHROMOSOME.getName());
		assertEquals("prefRef", GenomesOnly.PREF_REF.getName());
		assertEquals("exactly four values", 4, GenomesOnly.values().length);
	}

	/** A value is looked up by its configuration name, case insensitively, and an unknown one is null. */
	@Test
	public void testByName() {
		for (GenomesOnly value : GenomesOnly.values()) {
			assertSame(value, GenomesOnly.byName(value.getName()));
			assertSame(value, GenomesOnly.byName(value.getName().toUpperCase()));
		}
		assertNull(GenomesOnly.byName("nearlyComplete"));
		assertNull(GenomesOnly.byName(""));
	}

	/** The old boolean key must be gone rather than left beside its replacement. */
	@Test
	public void testTheBooleanKeyIsGone() throws Exception {
		for (GSConfigKey key : GSConfigKey.values()) {
			assertFalse("refseq.completeGenomesOnly still exists as " + key,
					"refseq.completeGenomesOnly".equals(key.getName()));
		}
	}
}
