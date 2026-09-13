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
package org.metagene.genestrip.store;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.nio.charset.StandardCharsets;
import java.util.Properties;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.junit.Test;
import org.metagene.genestrip.GSProject;

/**
 * Checks how {@link Database#load(java.io.InputStream, boolean)} reacts to a database file that names a
 * class which does not exist in the runtime library anymore.
 * <p>
 * This happens whenever a store implementation is dropped: a database written before the removal (say,
 * with the former {@code KMerSortedArray}) still carries that class name in its serialized stream, and
 * deserializing it raises a {@link ClassNotFoundException}. That is the same situation as a version
 * mismatch and must be reported as such - not as a bare missing-class error, which says nothing about
 * what to do next.
 */
public class DatabaseUnknownClassTest {
	/** Serializable stand-in for a store class that a later library version no longer has. */
	private static class Absent implements Serializable {
		private static final long serialVersionUID = 1L;
	}

	@Test
	public void testUnknownClassIsReportedAsVersionMismatch() throws IOException {
		byte[] db = buildDatabaseFile("0.9", "genestrip-core");

		try {
			Database.load(new ByteArrayInputStream(db), false);
			fail("Expected the load of a database naming an unknown class to fail.");
		} catch (InvalidDatabaseClassException e) {
			// The missing class must be named, so that it is clear which part of the database is at fault.
			assertTrue(e.getMessage(), e.getMessage().contains("Absnet"));
			// The config info is read even though the failure happens earlier in the file - only that
			// makes the message below able to name the version that wrote the database.
			assertEquals("0.9", e.getConfigInfo().getProperty(GSProject.GENESTRIP_VERSION));
			String message = e.toRuntimeException().getMessage();
			assertTrue(message, message.contains("0.9"));
		} catch (ClassNotFoundException e) {
			fail("The unknown class was reported as a bare " + e.getClass().getSimpleName()
					+ ", which does not tell the user that the database has to be regenerated: " + e.getMessage());
		}
	}

	// Builds a database file whose db entry names a class that does not exist, by serializing 'Absent'
	// and renaming it in the resulting bytes (to a name of the same length, so that the stream stays
	// well-formed). The config entry comes last, just as Database.save() writes it.
	private byte[] buildDatabaseFile(String version, String title) throws IOException {
		ByteArrayOutputStream serialized = new ByteArrayOutputStream();
		try (ObjectOutputStream oos = new ObjectOutputStream(serialized)) {
			oos.writeObject(new Absent());
		}
		byte[] entry = rename(serialized.toByteArray(), "Absent", "Absnet");

		Properties configInfo = new Properties();
		configInfo.setProperty(GSProject.GENESTRIP_VERSION, version);
		configInfo.setProperty(GSProject.GENESTRIP_TITLE, title);

		ByteArrayOutputStream result = new ByteArrayOutputStream();
		try (ZipOutputStream zos = new ZipOutputStream(result)) {
			zos.putNextEntry(new ZipEntry(Database.DB_FILE));
			zos.write(entry);
			zos.closeEntry();
			zos.putNextEntry(new ZipEntry(Database.CONFIG_INFO_FILE));
			configInfo.store(zos, null);
			zos.closeEntry();
		}
		return result.toByteArray();
	}

	// Replaces the first occurrence of 'from' by 'to' (both of equal length) in the given bytes.
	private byte[] rename(byte[] bytes, String from, String to) {
		byte[] fromBytes = from.getBytes(StandardCharsets.US_ASCII);
		byte[] toBytes = to.getBytes(StandardCharsets.US_ASCII);
		assertEquals("The replacement must not change the length of the stream.", fromBytes.length, toBytes.length);

		for (int i = 0; i <= bytes.length - fromBytes.length; i++) {
			boolean found = true;
			for (int j = 0; j < fromBytes.length; j++) {
				if (bytes[i + j] != fromBytes[j]) {
					found = false;
					break;
				}
			}
			if (found) {
				System.arraycopy(toBytes, 0, bytes, i, toBytes.length);
				return bytes;
			}
		}
		fail("Could not find '" + from + "' in the serialized stream.");
		return bytes;
	}
}
