package org.metagene.genestrip.goals.refseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.util.StreamProvider;

public class CategoriesGoal extends ObjectGoal<Set<RefSeqCategory>[], GSProject> {
	@SafeVarargs
	public CategoriesGoal(GSProject project, String name, Goal<GSProject>... deps) {
		super(project, name, deps);
	}
	
	@Override
	public void makeThis() {
		try {
			Set<RefSeqCategory>[] res = readFromFile(getProject().getCategoriesFile());
			set(res);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected Set<RefSeqCategory>[] readFromFile(File file) throws IOException {
		@SuppressWarnings("unchecked")
		Set<RefSeqCategory>[] res = new Set[2];
		res[0] = new HashSet<RefSeqCategory>();
		res[1] = new HashSet<RefSeqCategory>();

		InputStream stream = StreamProvider.getInputStreamForFile(file);

		try (BufferedReader br = new BufferedReader(new InputStreamReader(stream, StandardCharsets.UTF_8))) {
			String line = null;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				// Ignore whole line comments.
				if (!line.startsWith("#")) {
					// Ignore end of line comments.
					int commentIndex = line.indexOf('#');
					if (commentIndex != -1) {
						line = line.substring(0, commentIndex);
					}
					String category = line.trim();
					boolean additional = category.startsWith("+");
					if (additional) {
						category = category.substring(1);
					}
					RefSeqCategory cat = RefSeqCategory.fromDiregoryString(category);
					res[additional ? 1 : 0].add(cat);
				}
			}
		}
		
		res[1].addAll(res[0]);
		
		res[0] = Collections.unmodifiableSet(res[0]);
		res[1] = Collections.unmodifiableSet(res[1]);
		
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Categories for store: " + res[0]);
			getLogger().info("Categories for ancestor update: " + res[1]);
		}

		return res;
	}
}
