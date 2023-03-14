package org.metagene.genestrip.util;

public class CGAT {
	public static final byte[] CGAT_TO_UPPER_CASE = new byte[256];
	{
		for (int i = 0; i < CGAT_TO_UPPER_CASE.length; i++) {
			CGAT_TO_UPPER_CASE[i] = (byte) i;
		}
		CGAT_TO_UPPER_CASE['c'] = 'C';
		CGAT_TO_UPPER_CASE['g'] = 'G';
		CGAT_TO_UPPER_CASE['a'] = 'A';
		CGAT_TO_UPPER_CASE['t'] = 'T';
	}

	public byte cgatToUpperCase(byte bite) {
		return CGAT_TO_UPPER_CASE[bite];
	}

	public static boolean isCGAT(byte c) {
		return c == 'C' || c == 'G' || c == 'A' || c == 'T';
	}
}
