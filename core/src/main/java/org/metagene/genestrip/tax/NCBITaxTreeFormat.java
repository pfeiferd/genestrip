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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import me.tongfei.progressbar.ProgressBar;
import org.metagene.genestrip.io.BufferedLineReader;
import org.metagene.genestrip.io.ByteCountingInputStream;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.tax.TaxTree.TaxIdNode;
import org.metagene.genestrip.util.ByteArrayUtil;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;

/**
 * The NCBI taxonomy dump: the {@code nodes.dmp} and {@code names.dmp} files, read straight from the
 * archive {@code taxdmp.zip} as NCBI publishes it, or else from the two files themselves.
 */
public class NCBITaxTreeFormat implements TaxTreeFormat {
	private static final int MAX_LINE_SIZE = 4096;

	/** File name of the NCBI taxonomy nodes dump. */
	public static final String NODES_DMP = "nodes.dmp";
	/** File name of the NCBI taxonomy names dump. */
	public static final String NAMES_DMP = "names.dmp";
	/** File name of the NCBI taxonomy dump archive holding both. */
	public static final String TAX_DMP_ZIP = "taxdmp.zip";

	/**
	 * Reads {@code taxdmp.zip} in the given directory if it is there, and else {@code nodes.dmp} and
	 * {@code names.dmp}. The archive comes first because files extracted from an earlier archive
	 * may be outdated.
	 */
	@Override
	public TaxIdNode read(File path, TaxTree.Builder builder, boolean progressBar) throws IOException {
		File zipFile = new File(path, TAX_DMP_ZIP);
		if (zipFile.exists()) {
			return readZip(zipFile, builder, progressBar);
		}
		TaxIdNode root;
		File nodesFile = new File(path, NODES_DMP);
		File namesFile = new File(path, NAMES_DMP);
		if (progressBar) {
			try (StreamingResource.StreamAccess sa = new StreamingFileResource(nodesFile, true).openStream()) {
				try (ProgressBar pb = GSProgressBarCreator.newGSProgressBar("taxnodes", sa, null, true)) {
					root = readNodesFromStream(sa.getInputStream(), builder);
				}
			}
			try (StreamingResource.StreamAccess sa = new StreamingFileResource(namesFile, true).openStream()) {
				try (ProgressBar pb = GSProgressBarCreator.newGSProgressBar("taxnames", sa, null, true)) {
					readNamesFromStream(sa.getInputStream(), builder);
				}
			}
		}
		else {
			try (InputStream is = StreamProvider.getInputStreamForFile(nodesFile)) {
				root = readNodesFromStream(is, builder);
			}
			try (InputStream is = StreamProvider.getInputStreamForFile(namesFile)) {
				readNamesFromStream(is, builder);
			}
		}
		return root;
	}

	/**
	 * Reads {@code nodes.dmp} and {@code names.dmp} from the archive without extracting them.
	 *
	 * @param zipFile     the archive
	 * @param builder     the builder of the tree being read
	 * @param progressBar whether to show a progress bar while reading the entries
	 * @return the root node
	 * @throws IOException if reading the archive fails or an entry is missing
	 */
	protected TaxIdNode readZip(File zipFile, TaxTree.Builder builder, boolean progressBar) throws IOException {
		try (ZipFile zip = new ZipFile(zipFile)) {
			TaxIdNode root;
			try (EntryAccess ea = new EntryAccess(zip, NODES_DMP)) {
				try (ProgressBar pb = progressBar ? GSProgressBarCreator.newGSProgressBar("taxnodes", ea, null, true) : null) {
					root = readNodesFromStream(ea.getInputStream(), builder);
				}
			}
			try (EntryAccess ea = new EntryAccess(zip, NAMES_DMP)) {
				try (ProgressBar pb = progressBar ? GSProgressBarCreator.newGSProgressBar("taxnames", ea, null, true) : null) {
					readNamesFromStream(ea.getInputStream(), builder);
				}
			}
			return root;
		}
	}

	/**
	 * Reads a {@code names.dmp} stream and assigns names to the already loaded nodes,
	 * preferring scientific names.
	 *
	 * @param stream  the {@code names.dmp} input stream to read
	 * @param builder the builder of the tree being read
	 * @throws IOException if reading the stream fails
	 */
	protected void readNamesFromStream(InputStream stream, TaxTree.Builder builder) throws IOException {
		try (BufferedLineReader br = new BufferedLineReader(stream)) {
			int size;
			byte[] target = new byte[MAX_LINE_SIZE];

			while ((size = br.nextLine(target)) > 0) {
				boolean scientific = ByteArrayUtil.indexOf(target, 0, size, "scientific name") != -1;
				int a = ByteArrayUtil.indexOf(target, 0, size, '|');
				int b = ByteArrayUtil.indexOf(target, a + 1, size, '|');
				if (a != -1 && b != -1) {
					TaxIdNode node = builder.getNode(target, 0, a - 1, false);
					if (node != null) {
						if (node.getName() == null || scientific) {
							builder.setName(node, new String(target, a + 2, b - a - 3, StandardCharsets.UTF_8));
						}
					}
				}
			}
		}
	}

	/**
	 * Reads a {@code nodes.dmp} stream, creating nodes, linking them to their parents
	 * and setting their ranks.
	 *
	 * @param stream  the {@code nodes.dmp} input stream to read
	 * @param builder the builder of the tree being read
	 * @return the root node (tax id "1")
	 * @throws IOException if reading the stream fails
	 */
	protected TaxIdNode readNodesFromStream(InputStream stream, TaxTree.Builder builder) throws IOException {
		TaxIdNode res = null;
		int size;
		byte[] target = new byte[MAX_LINE_SIZE];

		try (BufferedLineReader br = new BufferedLineReader(stream)) {
			while ((size = br.nextLine(target)) > 0) {
				int a = ByteArrayUtil.indexOf(target, 0, size, '|');
				int b = ByteArrayUtil.indexOf(target, a + 1, size, '|');
				int c = ByteArrayUtil.indexOf(target, b + 1, size, '|');
				if (a != -1 && b != -1) {
					TaxIdNode nodeA = builder.getNode(target, 0, a - 1, true);
					TaxIdNode nodeB = builder.getNode(target, a + 2, b - 1, true);

					if (nodeA != nodeB) {
						builder.addSubNode(nodeB, nodeA);
					}
					builder.setRank(nodeA, Rank.getRankFromBytes(target, b + 2, c - 1));

					if (nodeA == nodeB && "1".equals(nodeA.taxId)) {
						res = nodeA;
					}
				}
			}
		}

		return res;
	}

	/**
	 * An entry of the archive as a stream that counts what was read from it, for the progress bar.
	 */
	protected static class EntryAccess implements StreamingResource.StreamAccess {
		private final ByteCountingInputStream inputStream;
		private final long size;

		/**
		 * Opens the named entry of the archive.
		 *
		 * @param zip  the archive
		 * @param name the name of the entry
		 * @throws IOException if the entry is missing or cannot be opened
		 */
		public EntryAccess(ZipFile zip, String name) throws IOException {
			ZipEntry entry = zip.getEntry(name);
			if (entry == null) {
				throw new IOException("Missing entry " + name + " in " + zip.getName());
			}
			size = entry.getSize();
			inputStream = new ByteCountingInputStream(zip.getInputStream(entry));
		}

		@Override
		public InputStream getInputStream() {
			return inputStream;
		}

		@Override
		public long getBytesRead() {
			return inputStream.getBytesRead();
		}

		@Override
		public long getSize() {
			return size;
		}
	}
}
