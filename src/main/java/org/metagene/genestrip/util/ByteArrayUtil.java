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
package org.metagene.genestrip.util;

import java.io.PrintStream;

public class ByteArrayUtil {
	public static int indexOf(byte[] array, int start, int end, char c) {
		for (int i = start; i < end; i++) {
			if (array[i] == c) {
				return i;
			}
		}
		return -1;
	}

	public static boolean equals(byte[] array, int start, int end, String str) {
		if (str.length() != end - start) {
			return false;
		}
		for (int i = start, j = 0; i < end; i++,j++) {
			if (array[i] != str.charAt(j)) {
				return false;
			}
		}
		return true;
	}

	public static int indexOf(byte[] array, int start, int end, String str) {
		int len = str.length();
		for (int i = start; i <= end - len; i++) {
			boolean found = true;
			for (int j = 0; j < len; j++) {
				if (array[i + j] != str.charAt(j)) {
					found = false;
					break;
				}
			}
			if (found) {
				return i;
			}
		}
		return -1;
	}

	public static boolean startsWith(byte[] array, int start, String str) {
		int len = str.length();
		if (array.length < start + len) {
			return false;
		}
		for (int j = 0; j < len; j++) {
			if (array[start + j] != str.charAt(j)) {
				return false;
			}
		}
		return true;
	}

	public static String toString(byte[] array) {
		int i;
		for (i = 0; i < array.length; i++) {
			if (array[i] == 0) {
				break;
			}
		}
		return new String(array, 0, i);
	}

	public static void println(byte[] array, PrintStream ps) {
		println(array, 0, array.length, ps);
	}

	public static void print(final byte[] array, final PrintStream ps) {
		print(array, 0, array.length, ps);
	}

	public static void println(byte[] array, int start, int end, PrintStream ps) {
		print(array, start, end, ps);
		ps.println();
	}

	public static void print(byte[] array, int start, int end, PrintStream ps) {
		for (int i = start; i < array.length && i < end; i++) {
			if (array[i] == 0) {
				break;
			}
			ps.print((char) array[i]);
		}
	}

	public static boolean checkDigits(byte[] data, int start, int end) {
		for (int i = start; i <= end; i++) {
			if (data[i] < '0' || data[i] > '9') {
				return false;
			}
		}
		return true;
	}

	public static int byteArrayToInt(byte[] data, int start, int end) {
		int result = 0;
		for (int i = start; i < end; i++) {
			int digit = data[i] - '0';
			if ((digit < 0) || (digit > 9)) {
				throw new NumberFormatException("For string byte array / string " + new String(data));
			}
			result *= 10;
			result += digit;
		}
		return result;
	}

	public static int intToByteArray(int value, final byte[] data, int pos) {
		if (value == 0) {
			data[pos++] = '0';
			return pos;
		}

		if (value < 0) {
			data[pos++] = '-';
		}
		int start = pos;
		while (value != 0) {
			data[pos++] = (byte) ('0' + (value % 10));
			value /= 10;
		}
		// Reorder:
		int end = pos - 1;
		byte h;
		while (start < end) {
			h = data[start];
			data[start] = data[end];
			data[end] = h;
			start++;
			end--;
		}
		return pos;
	}
}
