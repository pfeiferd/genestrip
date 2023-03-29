package org.metagene.genestrip.util;

public class ByteArrayToString {
	public static String toString(byte[] array) {
		int i;
		for (i = 0; i < array.length; i++) {
			if (array[i] == 0) {
				break;
			}
		}
		return new String(array, 0, i);
	}
}
