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
package org.metagene.genestrip.make;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.net.InetSocketAddress;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.codec.digest.DigestUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpServer;

/**
 * Tests that an HTTP download cut short is noticed and resumed, with and without a length announced by
 * the server, and that nothing appears under the final name before the download is complete.
 */
public class FileDownloadGoalTest {
	@Rule
	public final TemporaryFolder folder = new TemporaryFolder();

	private HttpServer server;
	private byte[] payload;
	/** How many responses are cut short, at half their length. */
	private int cuts;
	/** Whether the length is announced; if not, the response is chunked. */
	private boolean announceLength;
	/** Whether range requests are honoured. */
	private boolean ranges;
	/** Whether the remote file is replaced by a new version after the first response. */
	private boolean changeAfterFirst;
	/** The version of the remote file, sent as its ETag. */
	private String etag;
	private int requests;
	/** The Range header of every request, in order. */
	private final List<String> rangesAsked = new ArrayList<>();

	@Before
	public void setUp() throws IOException {
		Random random = new Random(42);
		ByteArrayOutputStream bytes = new ByteArrayOutputStream();
		try (GZIPOutputStream gzip = new GZIPOutputStream(bytes)) {
			byte[] text = new byte[1 << 20];
			for (int i = 0; i < text.length; i++) {
				text[i] = (byte) "ACGT".charAt(random.nextInt(4));
			}
			gzip.write(text);
		}
		payload = bytes.toByteArray();
		etag = "\"v1\"";
		announceLength = true;
		ranges = true;
		server = HttpServer.create(new InetSocketAddress("127.0.0.1", 0), 0);
		server.createContext("/", this::handle);
		server.start();
	}

	@After
	public void tearDown() {
		server.stop(0);
	}

	private void handle(HttpExchange exchange) throws IOException {
		requests++;
		String path = exchange.getRequestURI().getPath();
		if (!path.endsWith("/data.gz") && !path.endsWith("/data.txt")) {
			exchange.sendResponseHeaders(404, -1);
			exchange.close();
			return;
		}
		int start = 0;
		String range = exchange.getRequestHeaders().getFirst("Range");
		rangesAsked.add(range);
		String ifRange = exchange.getRequestHeaders().getFirst("If-Range");
		if (ifRange != null && !ifRange.equals(etag)) {
			// The version asked for is gone: the whole of the current one is sent instead.
			range = null;
		}
		if (ranges && range != null) {
			start = Integer.parseInt(range.substring("bytes=".length(), range.indexOf('-')));
		}
		exchange.getResponseHeaders().set("ETag", etag);
		int length = payload.length - start;
		int send = cuts > 0 ? length / 2 : length;
		boolean cut = cuts > 0;
		cuts--;
		if (ranges && range != null) {
			exchange.getResponseHeaders().set("Content-Range", "bytes " + start + "-" + (payload.length - 1) + "/"
					+ (announceLength ? String.valueOf(payload.length) : "*"));
		}
		exchange.sendResponseHeaders(ranges && range != null ? 206 : 200, announceLength ? length : 0);
		OutputStream out = exchange.getResponseBody();
		try {
			out.write(payload, start, send);
			if (cut && announceLength) {
				// Fewer bytes than announced: the connection is dropped.
				out.flush();
				throw new IOException("cut");
			}
			out.close();
		} catch (IOException e) {
			exchange.close();
		}
		if (changeAfterFirst) {
			changeAfterFirst = false;
			payload = Arrays.copyOf(payload, payload.length);
			payload[payload.length / 4] ^= 1;
			etag = "\"v2\"";
		}
	}

	private Goal download(String name, String md5) {
		File file = new File(folder.getRoot(), name);
		return new Goal(new TestProject(), file, md5, "http://127.0.0.1:" + server.getAddress().getPort());
	}

	@Test
	public void testCompleteDownload() throws IOException {
		Goal goal = download("data.gz", null);
		goal.make();
		assertArrayEquals(payload, Files.readAllBytes(goal.getFile().toPath()));
		assertEquals(1, requests);
	}

	@Test
	public void testCutWithLengthIsResumed() throws IOException {
		cuts = 2;
		Goal goal = download("data.gz", null);
		goal.make();
		assertArrayEquals(payload, Files.readAllBytes(goal.getFile().toPath()));
		assertEquals(3, requests);
		assertFalse(goal.partFile(goal.getFile()).exists());
	}

	@Test
	public void testCutWithoutLengthIsResumed() throws IOException {
		// A chunked response that ends cleanly but early: only the gzip stream shows it.
		announceLength = false;
		cuts = 2;
		Goal goal = download("data.gz", null);
		goal.make();
		assertArrayEquals(payload, Files.readAllBytes(goal.getFile().toPath()));
		assertEquals(3, requests);
	}

	@Test
	public void testServerIgnoringRangesStartsOver() throws IOException {
		ranges = false;
		cuts = 1;
		Goal goal = download("data.gz", null);
		goal.make();
		assertArrayEquals(payload, Files.readAllBytes(goal.getFile().toPath()));
		assertEquals(2, requests);
	}

	@Test
	public void testPartOfEarlierProcessIsNotResumed() throws IOException {
		// Nothing says which version of the remote file a leftover part belongs to.
		Goal goal = download("data.gz", null);
		Files.write(goal.partFile(goal.getFile()).toPath(), new byte[1000]);
		goal.make();
		assertArrayEquals(payload, Files.readAllBytes(goal.getFile().toPath()));
		assertEquals("bytes=0-", rangesAsked.get(0));
	}

	@Test
	public void testChangedRemoteFileStartsOver() throws IOException {
		byte[] first = payload;
		changeAfterFirst = true;
		cuts = 1;
		Goal goal = download("data.gz", null);
		goal.make();
		// The remainder of the new version must not be appended to the beginning of the old one.
		assertFalse(Arrays.equals(first, payload));
		assertArrayEquals(payload, Files.readAllBytes(goal.getFile().toPath()));
		assertEquals(2, requests);
	}

	@Test
	public void testNonGzipWithoutLengthIsTakenAsItEnds() throws IOException {
		// Without a length and without gzip there is nothing to check a clean end against.
		announceLength = false;
		Goal goal = download("data.txt", null);
		goal.make();
		assertArrayEquals(payload, Files.readAllBytes(goal.getFile().toPath()));
		assertEquals(1, requests);
	}

	@Test
	public void testCheckSumIsVerified() throws IOException {
		cuts = 1;
		Goal goal = download("data.gz", DigestUtils.md5Hex(payload));
		goal.make();
		assertArrayEquals(payload, Files.readAllBytes(goal.getFile().toPath()));
		assertTrue(goal.isMade());
	}

	@Test
	public void testCutOnEveryTryFails() {
		cuts = Integer.MAX_VALUE;
		Goal goal = download("data.gz", null);
		try {
			goal.make();
			fail("an incomplete download was accepted");
		} catch (RuntimeException e) {
			// Expected.
		}
		assertFalse(goal.getFile().exists());
		assertFalse(goal.partFile(goal.getFile()).exists());
	}

	@Test
	public void testMissingFile() {
		Goal goal = download("missing.gz", null);
		try {
			goal.make();
			fail("a missing file was accepted");
		} catch (RuntimeException e) {
			assertTrue(String.valueOf(e.getCause()), e.getCause() instanceof FileNotFoundException);
		}
	}

	@Test
	public void testContentRangeParsing() {
		assertEquals(100, FileDownloadGoal.rangeStart("bytes 100-199/1000"));
		assertEquals(1000, FileDownloadGoal.rangeTotal("bytes 100-199/1000"));
		assertEquals(-1, FileDownloadGoal.rangeTotal("bytes 100-199/*"));
		assertEquals(-1, FileDownloadGoal.rangeStart(null));
	}

	private static final class TestProject extends Project {
		private TestProject() {
			super("test");
		}

		@Override
		protected ConfigKey[] getConfigKeys() {
			return new ConfigKey[] { FileDownloadGoal.CHECK_SUM_CACHE_FILE };
		}
	}

	private static final class Goal extends FileDownloadGoal<TestProject> {
		private final File file;
		private final String md5;
		private final String baseURL;

		private Goal(TestProject project, File file, String md5, String baseURL) {
			super(project, new GoalKey.DefaultGoalKey("download"));
			this.file = file;
			this.md5 = md5;
			this.baseURL = baseURL;
		}

		@Override
		public List<File> getFiles() {
			return Collections.singletonList(file);
		}

		@Override
		protected String getMD5CheckSum(File file) {
			return md5;
		}

		@Override
		protected String getHttpBaseURL() {
			return baseURL;
		}

		@Override
		protected String getFTPBaseURL() {
			return null;
		}

		@Override
		protected boolean isUseHttp() {
			return true;
		}

		@Override
		protected boolean isIgnoreMissingFiles() {
			return false;
		}

		@Override
		protected int getMaxDownloadTries() {
			return 3;
		}

		@Override
		protected String getFTPDir(File file) {
			return "/dir";
		}
	}
}
