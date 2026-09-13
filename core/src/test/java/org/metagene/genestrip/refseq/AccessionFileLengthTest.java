/*
 * Genestrip
 */
package org.metagene.genestrip.refseq;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import org.junit.Test;
import org.metagene.genestrip.GSConfigKey.RefSeqStatus;
import org.metagene.genestrip.GSConfigKey.SeqType;
import org.metagene.genestrip.io.StreamingFileResource;

/**
 * Tests that the length column of the accession catalog reaches a subclass intact.
 * <p>
 * The length is the last column of the line, and the line reader writes the terminator into the
 * buffer and counts it, so handing the line's end on as the column's end gives a subclass
 * {@code "7768\n"} instead of {@code "7768"}. That parses as nothing, which cost every row of a
 * three gigabyte catalog before it was caught - the pass recognised no assembly at all, and nothing
 * in it failed loudly enough to say why.
 */
public class AccessionFileLengthTest {

	private static final String CATALOG = "9\tBuchnera aphidicola\tNC_001911.1\tbacteria|complete|plasmid\tna\t7768\n"
			+ "9\tBuchnera aphidicola\tNC_004843.1\tbacteria|complete|plasmid\tna\t2308\n"
			+ "2\tBacteria\tWP_000002109.1\tbacteria|complete\tUNKNOWN\t93\n"
			+ "562\tEscherichia coli\tNC_000913.3\tbacteria|complete\tna\t4641652\n";

	/** A catalog whose last line carries no terminator at all must parse like any other. */
	private static final String NO_TRAILING_NEWLINE = CATALOG.substring(0, CATALOG.length() - 1);

	private List<long[]> lengthsOf(String catalog) throws IOException {
		File f = File.createTempFile("catalog", ".gz");
		f.deleteOnExit();
		try (GZIPOutputStream out = new GZIPOutputStream(new FileOutputStream(f))) {
			out.write(catalog.getBytes("UTF-8"));
		}
		final List<long[]> seen = new ArrayList<long[]>();
		AccessionFileProcessor processor = new AccessionFileProcessor(
				new java.util.HashSet<RefSeqCategory>(Arrays.asList(RefSeqCategory.BACTERIA)), SeqType.GENOMIC,
				Arrays.asList(RefSeqStatus.values())) {
			@Override
			protected void handleEntry(byte[] target, int taxIdEnd, int accessionStart, int accessionEnd) {
				// Not called: the overload below is.
			}

			@Override
			protected void handleEntry(byte[] target, int taxIdEnd, int accessionStart, int accessionEnd,
					int lengthStart, int lengthEnd) {
				seen.add(new long[] { parseLength(target, 0, taxIdEnd),
						parseLength(target, lengthStart, lengthEnd) });
			}

			@Override
			protected boolean containsCategory(byte[] target, int start, int end, RefSeqCategory[] categories) {
				// The catalog's categories are not what is under test here.
				return true;
			}
		};
		processor.processCatalog(new StreamingFileResource(f));
		return seen;
	}

	/** Every genomic row's length must arrive as the number it is, not as the number plus a newline. */
	@Test
	public void testLengthsParse() throws IOException {
		List<long[]> seen = lengthsOf(CATALOG);
		assertEquals("three genomic rows, the protein one filtered", 3, seen.size());
		assertEquals(9, seen.get(0)[0]);
		assertEquals("the length must parse, not come back as -1", 7768, seen.get(0)[1]);
		assertEquals(2308, seen.get(1)[1]);
		assertEquals(562, seen.get(2)[0]);
		assertEquals(4641652, seen.get(2)[1]);
		for (long[] row : seen) {
			assertTrue("no row may fail to parse", row[1] > 0);
		}
	}

	/** A file whose final line has no terminator must behave the same. */
	@Test
	public void testLastLineWithoutTerminator() throws IOException {
		List<long[]> seen = lengthsOf(NO_TRAILING_NEWLINE);
		assertEquals(3, seen.size());
		assertEquals(4641652, seen.get(2)[1]);
	}
}
