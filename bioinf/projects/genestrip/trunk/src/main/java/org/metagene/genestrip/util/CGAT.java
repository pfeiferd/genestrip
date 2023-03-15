package org.metagene.genestrip.util;

public class CGAT {
	private static final byte[] CGAT_TO_UPPER_CASE = new byte[Byte.MAX_VALUE + 1];
	{
		for (int i = 0; i < CGAT_TO_UPPER_CASE.length; i++) {
			CGAT_TO_UPPER_CASE[i] = (byte) i;
		}
		CGAT_TO_UPPER_CASE['c'] = 'C';
		CGAT_TO_UPPER_CASE['g'] = 'G';
		CGAT_TO_UPPER_CASE['a'] = 'A';
		CGAT_TO_UPPER_CASE['t'] = 'T';
	}

	public static byte cgatToUpperCase(byte c) {
		if (c == 'C' || c == 'G' || c == 'A' || c == 'T') {
			return CGAT_TO_UPPER_CASE[c];
		}
		return c;
	}

	public static boolean isCGAT(byte c) {
		return c == 'C' || c == 'G' || c == 'A' || c == 'T';
	}
}
