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
package org.metagene.genestrip.accuracy;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.fasta.AbstractFastaReader;
import org.metagene.genestrip.io.BufferedLineReader;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.util.ByteArrayUtil;

public class AccessionNumber2TaxidGoal extends ObjectGoal<Map<String, String>, GSProject> {
	public static final String ACCESSION_MAP_FILE = "nucl_gb.accession2taxid.gz";
	public static final String OLD_ACCESSION_MAP_FILE = "dead_nucl.accession2taxid.gz";

	private static final String MARKER = "|ref|";

	private final File fastaFile;

	@SafeVarargs
	public AccessionNumber2TaxidGoal(GSProject project, String name, File fastaFile, Goal<GSProject>... dependencies) {
		super(project, name, dependencies);
		this.fastaFile = fastaFile;
	}

	@Override
	public void makeThis() {
		final Map<String, String> accesion2TaxidMap = new HashMap<String, String>();
		Set<String> accNumbers = new HashSet<String>();
		AbstractFastaReader fastaReader1 = new AbstractFastaReader(getProject().getConfig().getMaxReadSizeBytes()) {
			@Override
			protected void infoLine() throws IOException {
				String acc = getAccessionNumberFromInfoLine(target, size);
				if (acc != null) {
					accNumbers.add(acc);
				}
			}

			@Override
			protected void dataLine() throws IOException {
			}
		};
		try {
			fastaReader1.readFasta(fastaFile);
			fillAccessionNumbersToTaxIdsMap(ACCESSION_MAP_FILE, accesion2TaxidMap, accNumbers);
			fillAccessionNumbersToTaxIdsMap(OLD_ACCESSION_MAP_FILE, accesion2TaxidMap, accNumbers);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		set(accesion2TaxidMap);
	}

	protected void fillAccessionNumbersToTaxIdsMap(String fileName, Map<String, String> accesion2TaxidMap,
			Set<String> accNumbers) throws IOException {

		InputStream in = StreamProvider
				.getInputStreamForFile(new File(getProject().getConfig().getCommonDir(), fileName));

		BufferedLineReader br = new BufferedLineReader(in);
		int size;
		byte[] target = new byte[2048];
		br.nextLine(target);
		while ((size = br.nextLine(target)) > 0) {
			int tab1 = ByteArrayUtil.indexOf(target, 0, size, '\t');
			int tab2 = ByteArrayUtil.indexOf(target, tab1 + 1, size, '\t');
			String acc = new String(target, tab1 + 1, tab2 - tab1 - 1);
			if (accNumbers.contains(acc)) {
				int tab3 = ByteArrayUtil.indexOf(target, tab2 + 1, size, '\t');
				String tax = new String(target, tab2 + 1, tab3 - tab2 - 1);
				accesion2TaxidMap.put(acc, tax);
			}
		}
	}

	public static String getAccessionNumberFromInfoLine(byte[] target, int size) {
		int pos1 = ByteArrayUtil.indexOf(target, 1, size, MARKER);
		if (pos1 != -1) {
			pos1 += MARKER.length();
			int pos2 = ByteArrayUtil.indexOf(target, pos1, size, '|');
			if (pos2 != -1) {
				return new String(target, pos1, pos2 - pos1).trim();
			}
		}
		return null;
	}
}
