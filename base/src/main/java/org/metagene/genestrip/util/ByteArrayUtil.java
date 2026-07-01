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
import java.nio.charset.StandardCharsets;

import java.io.PrintStream;

/**
 * Static helpers for working with byte arrays that hold ASCII / UTF-8 text: searching, comparing,
 * printing, and converting between integers and their decimal ASCII representation.
 */
public class ByteArrayUtil {
	// Static utility class - not meant to be instantiated.
	private ByteArrayUtil() {
	}

	/**
	 * Returns the index of the first occurrence of {@code c} within {@code [start, end)}, or
	 * {@code -1} if it is not found.
	 *
	 * @param array the byte array to search
	 * @param start the index (inclusive) at which to start searching
	 * @param end   the index (exclusive) at which to stop searching
	 * @param c     the character to look for
	 * @return the index of the first matching byte, or {@code -1} if none is found
	 */
	public static int indexOf(byte[] array, int start, int end, char c) {
		for (int i = start; i < end; i++) {
			if (array[i] == c) {
				return i;
			}
		}
		return -1;
	}

	/**
	 * Returns whether the bytes in {@code [start, end)} equal, character by character, the given
	 * string.
	 *
	 * @param array the byte array to compare
	 * @param start the index (inclusive) at which the comparison starts
	 * @param end   the index (exclusive) at which the comparison stops
	 * @param str   the string to compare against
	 * @return {@code true} if the byte range equals the string, {@code false} otherwise
	 */
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

	/**
	 * Returns the start index of the first occurrence of {@code str} within {@code [start, end)}, or
	 * {@code -1} if it is not found.
	 *
	 * @param array the byte array to search
	 * @param start the index (inclusive) at which to start searching
	 * @param end   the index (exclusive) at which to stop searching
	 * @param str   the string to look for
	 * @return the start index of the first match, or {@code -1} if none is found
	 */
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

	/**
	 * Returns whether the array, starting at index {@code start}, begins with the given string.
	 *
	 * @param array the byte array to test
	 * @param start the index at which to start the comparison
	 * @param str   the prefix string to test for
	 * @return {@code true} if the array starts with the string at {@code start}, {@code false} otherwise
	 */
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

	/**
	 * Returns the UTF-8 string formed by the bytes up to (but excluding) the first null (0) byte, or
	 * the whole array if it contains no null byte.
	 *
	 * @param array the byte array to convert
	 * @return the decoded UTF-8 string
	 */
	public static String toString(byte[] array) {
		int i;
		for (i = 0; i < array.length; i++) {
			if (array[i] == 0) {
				break;
			}
		}
		return new String(array, 0, i, StandardCharsets.UTF_8);
	}

	/**
	 * Prints the whole array (up to the first null byte) followed by a line separator.
	 *
	 * @param array the byte array to print
	 * @param ps    the print stream to write to
	 */
	public static void println(byte[] array, PrintStream ps) {
		println(array, 0, array.length, ps);
	}

	/**
	 * Prints the whole array's characters up to the first null byte.
	 *
	 * @param array the byte array to print
	 * @param ps    the print stream to write to
	 */
	public static void print(final byte[] array, final PrintStream ps) {
		print(array, 0, array.length, ps);
	}

	/**
	 * Prints the characters in {@code [start, end)} followed by a line separator.
	 *
	 * @param array the byte array to print
	 * @param start the index (inclusive) at which to start printing
	 * @param end   the index (exclusive) at which to stop printing
	 * @param ps    the print stream to write to
	 */
	public static void println(byte[] array, int start, int end, PrintStream ps) {
		print(array, start, end, ps);
		ps.println();
	}

	/**
	 * Prints the characters in {@code [start, end)}, stopping early at the first null byte.
	 *
	 * @param array the byte array to print
	 * @param start the index (inclusive) at which to start printing
	 * @param end   the index (exclusive) at which to stop printing
	 * @param ps    the print stream to write to
	 */
	public static void print(byte[] array, int start, int end, PrintStream ps) {
		for (int i = start; i < array.length && i < end; i++) {
			if (array[i] == 0) {
				break;
			}
			ps.print((char) array[i]);
		}
	}

	/**
	 * Returns whether every byte in {@code [start, end)} is an ASCII decimal digit.
	 *
	 * @param data  the byte array to inspect
	 * @param start the index (inclusive) at which to start checking
	 * @param end   the index (exclusive) at which to stop checking
	 * @return {@code true} if all bytes in the range are ASCII digits, {@code false} otherwise
	 */
	public static boolean checkDigits(byte[] data, int start, int end) {
		for (int i = start; i < end; i++) {
			if (data[i] < '0' || data[i] > '9') {
				return false;
			}
		}
		return true;
	}

	/**
	 * Parses the ASCII decimal digits in {@code [start, end)} into an {@code int}.
	 *
	 * @param data  the byte array containing the digits
	 * @param start the index (inclusive) at which to start parsing
	 * @param end   the index (exclusive) at which to stop parsing
	 * @return the parsed integer value
	 * @throws NumberFormatException if a non-digit byte is encountered
	 */
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

	/**
	 * Writes the decimal ASCII representation of {@code value} (including a leading {@code '-'} for
	 * negative values) into {@code data} starting at {@code pos}.
	 *
	 * @param value the integer value to write
	 * @param data  the byte array to write into
	 * @param pos   the index at which to start writing
	 * @return the index in {@code data} just after the last written byte
	 */
	public static int intToByteArray(int value, final byte[] data, int pos) {
		if (value == 0) {
			data[pos++] = '0';
			return pos;
		}
		if (value < 0) {
			data[pos++] = '-';
		} else {
			// Extract digits in the negative range so Integer.MIN_VALUE (which has no positive
			// counterpart) is handled correctly.
			value = -value;
		}
		int start = pos;
		while (value != 0) {
			// value is <= 0 here, so value % 10 is in [-9, 0] and '0' - (value % 10) is the digit.
			data[pos++] = (byte) ('0' - (value % 10));
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
