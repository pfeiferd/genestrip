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
package org.metagene.genestrip.tax;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.fail;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;

/**
 * Tests that {@link NCBITaxTreeFormat} reads the tree straight from {@code taxdmp.zip}, that the
 * archive wins over extracted files beside it, and that the extracted files alone still work.
 */
public class NCBITaxTreeFormatTest {
	@Rule
	public final TemporaryFolder folder = new TemporaryFolder();

	private static final String NODES = "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tgenus\t|\n3\t|\t2\t|\tspecies\t|\n";
	private static final String NAMES = "1\t|\troot\t|\t\t|\tscientific name\t|\n"
			+ "2\t|\tGenus synonym\t|\t\t|\tsynonym\t|\n2\t|\tGenus\t|\t\t|\tscientific name\t|\n"
			+ "3\t|\tSpecies\t|\t\t|\tscientific name\t|\n" + "99\t|\tNot in nodes\t|\t\t|\tscientific name\t|\n";

	private void writeZip(File dir, String nodes, String names) throws IOException {
		try (ZipOutputStream zos = new ZipOutputStream(
				new FileOutputStream(new File(dir, NCBITaxTreeFormat.TAX_DMP_ZIP)))) {
			for (String[] entry : new String[][] { { NCBITaxTreeFormat.NODES_DMP, nodes },
					{ NCBITaxTreeFormat.NAMES_DMP, names } }) {
				if (entry[1] != null) {
					zos.putNextEntry(new ZipEntry(entry[0]));
					zos.write(entry[1].getBytes(StandardCharsets.UTF_8));
					zos.closeEntry();
				}
			}
		}
	}

	private void writeDmp(File dir, String nodes, String names) throws IOException {
		Files.write(new File(dir, NCBITaxTreeFormat.NODES_DMP).toPath(), nodes.getBytes(StandardCharsets.UTF_8));
		Files.write(new File(dir, NCBITaxTreeFormat.NAMES_DMP).toPath(), names.getBytes(StandardCharsets.UTF_8));
	}

	private void checkTree(TaxTree tree) {
		TaxIdNode species = tree.getNodeByTaxId("3");
		assertEquals("Species", species.getName());
		assertEquals(Rank.SPECIES, species.getRank());
		assertSame(tree.getNodeByTaxId("2"), species.getParent());
		// The scientific name wins over a synonym listed before it.
		assertEquals("Genus", tree.getNodeByTaxId("2").getName());
		assertSame(tree.getNodeByTaxId("1"), species.getParent().getParent());
		// A name for a tax id that nodes.dmp lacks creates no node.
		assertNull(tree.getNodeByTaxId("99"));
	}

	@Test
	public void testReadsFromZip() throws IOException {
		File dir = folder.newFolder();
		writeZip(dir, NODES, NAMES);
		checkTree(new TaxTree(dir, false));
		checkTree(new TaxTree(dir, true));
	}

	@Test
	public void testReadsFromDmpFiles() throws IOException {
		File dir = folder.newFolder();
		writeDmp(dir, NODES, NAMES);
		checkTree(new TaxTree(dir, false, new NCBITaxTreeFormat()));
	}

	@Test
	public void testZipWinsOverExtractedFiles() throws IOException {
		File dir = folder.newFolder();
		writeZip(dir, NODES, NAMES);
		// Outdated files from an earlier archive, where the species does not exist yet.
		writeDmp(dir, "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tgenus\t|\n", "1\t|\troot\t|\t\t|\tscientific name\t|\n");
		checkTree(new TaxTree(dir, false));
	}

	@Test
	public void testMissingEntry() throws IOException {
		File dir = folder.newFolder();
		writeZip(dir, NODES, null);
		try {
			new TaxTree(dir, false);
			fail("a missing names.dmp went unnoticed");
		} catch (RuntimeException e) {
			assertEquals(IOException.class, e.getCause().getClass());
		}
	}
}
