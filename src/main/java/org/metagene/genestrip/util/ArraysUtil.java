package org.metagene.genestrip.util;

import java.util.Arrays;

public class ArraysUtil {
	@SafeVarargs
	public static <T> T[] append(T[] array, T... values) {
	     T[] result = Arrays.copyOf(array, array.length + values.length);
	     for (int i = array.length; i < result.length; i++) {
	    	 result[i] = values[i - array.length];	    	 
	     }
	     return result;
	}
}
