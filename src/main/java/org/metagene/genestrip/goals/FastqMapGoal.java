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
package org.metagene.genestrip.goals;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingFileResource;
import org.metagene.genestrip.io.StreamingResource;
import org.metagene.genestrip.io.StreamingURLResource;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;

public class FastqMapGoal extends ObjectGoal<Map<String, List<StreamingResource>>, GSProject> {
	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#')
			.setDelimiter(' ').setRecordSeparator('\n').build();

	@SafeVarargs
	public FastqMapGoal(GSProject project, Goal<GSProject>... deps) {
		super(project, GSGoalKey.FASTQ_MAP, deps);
	}

	@Override
	protected void doMakeThis() {
		Map<String, List<StreamingResource>> map = createFastqMap(getProject().getKey(),
				getProject().getFastqResources(), getProject().getFastqMapFile());
		set(map);
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Derived fastq map: " + map);
		}
	}

	protected Map<String, List<StreamingResource>> createFastqMap(String key, String[] fastqs, String mapFilePath) {
		// Linked hash map preserve order of keys as entered.
		Map<String, List<StreamingResource>> resMap = new LinkedHashMap<String, List<StreamingResource>>();

		if (fastqs != null && fastqs.length > 0) {
			List<StreamingResource> resources = new ArrayList<StreamingResource>();
			for (String pathOrURL : fastqs) {
				List<StreamingResource> res = getResources(pathOrURL);
				if (res != null) {
					resources.addAll(res);
				} else if (getLogger().isWarnEnabled()) {
					getLogger().warn("Missing fastq resource: " + pathOrURL);
				}
			}
			if (!resources.isEmpty()) {
				if (key == null) {
					key = getProject().getFileBaseName(resources.get(0).getName());
				}
				resMap.put(key, resources);
			}
		}
		if (mapFilePath != null) {
			File csvFile = getProject().fastqFileFromPath(mapFilePath);
			if (csvFile != null) {
				try (CSVParser parser = FORMAT
						.parse(new InputStreamReader(StreamProvider.getInputStreamForFile(csvFile)))) {
					for (CSVRecord record : parser) {
						key = record.get(0);
						List<StreamingResource> resources = getResources(record.get(1));
						if (resources != null) {
							List<StreamingResource> l = resMap.get(key);
							if (l == null) {
								l = new ArrayList<StreamingResource>();
								resMap.put(key, l);
							}
							l.addAll(resources);
						} else if (getLogger().isWarnEnabled()) {
							getLogger().warn("Missing fastq resource: " + record.get(1));
						}
					}
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		}
		return resMap;
	}

	protected List<StreamingResource> getResources(String pathOrURL) {
		try {
			return Collections.singletonList(new StreamingURLResource(new URL(pathOrURL)));
		} catch (MalformedURLException e) {
			List<File> fastqs = getProject().fastqFilesFromPath(pathOrURL);
			if (fastqs != null) {
				List<StreamingResource> res = new ArrayList<StreamingResource>();
				for (File file : fastqs) {
					res.add(new StreamingFileResource(file));
				}
				return res;
			}
			return null;
		}
	}

}
