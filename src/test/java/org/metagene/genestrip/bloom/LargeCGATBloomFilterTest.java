package org.metagene.genestrip.bloom;

public class LargeCGATBloomFilterTest extends CGATBloomFilterTest {
	@Override
	protected boolean isTestLarge() {
		return true;
	}
}
