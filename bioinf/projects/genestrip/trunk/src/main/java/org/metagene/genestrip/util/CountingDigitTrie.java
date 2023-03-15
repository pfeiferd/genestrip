package org.metagene.genestrip.util;

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Map;

public class CountingDigitTrie {
	private CountingDigitTrie[] children;
	private long value;
	private String digits;

	public CountingDigitTrie() {
	}

	public String inc(String digits) {
		int index;
		int end = digits.length();
		CountingDigitTrie node = this, child;
		for (int i = 0; i < end; i++) {
			index = digits.charAt(i) - '0';
			if (node.children == null) {
				node.children = new CountingDigitTrie[10];
			}
			child = node.children[index];
			if (child == null) {
				child = new CountingDigitTrie();
				node.children[index] = child;
			}
		}
		node.value++;
		if (node.value == 1) {
			node.digits = digits;
		}
		
		return node.digits;
	}
	
	public String inc(byte[] seq, int start, int end) {
		int index;
		CountingDigitTrie node = this, child;
		for (int i = start; i <= end; i++, node = child) {
			index = seq[i] - '0';
			if (node.children == null) {
				node.children = new CountingDigitTrie[10];
			}
			child = node.children[index];
			if (child == null) {
				child = new CountingDigitTrie();
				node.children[index] = child;
			}
		}
		node.value++;
		if (node.value == 1) {
			node.digits = new String(seq, start, end - start + 1);
		}
		return node.digits;
	}

	public void collect(Map<String, Long> map) {
		if (value > 0) {
			map.put(digits, value);
		}
		if (children != null) {
			for (int k = 0; k < children.length; k++) {
				if (children[k] != null) {
					children[k].collect(map);
				}
			}
		}
	}
	
	public static void print(Map<String, Long> map, OutputStream out) {
		PrintStream pOut = new PrintStream(out);

		for (String taxid : map.keySet()) {
			pOut.print(taxid);
			pOut.print(';');
			pOut.println(map.get(taxid));
		}
	}	
}

