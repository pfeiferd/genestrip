package org.metagene.genestrip.util;

public class CGAT {
	private static final byte[] CGAT_TO_UPPER_CASE = new byte[256];

	static {
		for (int i = 0; i < CGAT_TO_UPPER_CASE.length; i++) {
			CGAT_TO_UPPER_CASE[i] = (byte) (i - 128);
		}
		CGAT_TO_UPPER_CASE[128 + 'c'] = 'C';
		CGAT_TO_UPPER_CASE[128 + 'g'] = 'G';
		CGAT_TO_UPPER_CASE[128 + 'a'] = 'A';
		CGAT_TO_UPPER_CASE[128 + 't'] = 'T';
	}

	public static byte cgatToUpperCase(byte c) {
		return CGAT_TO_UPPER_CASE[128 + c];
	}

	public static boolean isCGAT(byte c) {
		return c == 'C' || c == 'G' || c == 'A' || c == 'T';
	}
}
