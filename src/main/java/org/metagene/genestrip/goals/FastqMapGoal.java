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
import org.metagene.genestrip.io.StreamingResourceListStream;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.io.StreamingURLResource;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;

public class FastqMapGoal extends ObjectGoal<Map<String, StreamingResourceStream>, GSProject> {
	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#')
			.setDelimiter(' ').setRecordSeparator('\n').build();

	private boolean fastqType;

	@SafeVarargs
	public FastqMapGoal(GSProject project, boolean fastqType, Goal<GSProject>... deps) {
		super(project, fastqType ? GSGoalKey.FASTQ_MAP : GSGoalKey.FASTA_MAP, deps);
		this.fastqType = fastqType;
	}

	@Override
	protected void doMakeThis() {
		Map<String, StreamingResourceStream> map = createFastqMap(getProject().getKey(),
				getProject().getFastqResources(), getProject().getExtraResourcesKey(), getProject().getExtraResources(),
				getProject().getFastqMapFile());
		set(map);
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Derived fastq / fasta map: " + map);
		}
	}

	protected Map<String, StreamingResourceStream> createFastqMap(String key, String[] fastqs, String orKey,
			StreamingResourceStream otherResources, String mapFilePath) {
		// Linked hash map preserve order of keys as entered.
		Map<String, StreamingResourceStream> resMap = new LinkedHashMap<String, StreamingResourceStream>();

		if (fastqs != null && fastqs.length > 0) {
			StreamingResourceListStream resources = new StreamingResourceListStream();
			for (String pathOrURL : fastqs) {
				List<StreamingResource> res = getResources(pathOrURL);
				if (res != null) {
					resources.getList().addAll(res);
				} else if (getLogger().isWarnEnabled()) {
					getLogger().warn("Missing fastq or fasta resource: " + pathOrURL);
				}
			}
			if (!resources.getList().isEmpty()) {
				if (key == null) {
					key = getProject().getFileBaseName(resources.getList().get(0).getName());
				}
				resMap.put(key, resources);
			}
		}
		if (mapFilePath != null) {
			File csvFile = fastqType ? getProject().fastqFileFromPath(mapFilePath) : getProject().fastaFileFromPath(mapFilePath);
			if (csvFile != null) {
				try (CSVParser parser = FORMAT
						.parse(new InputStreamReader(StreamProvider.getInputStreamForFile(csvFile)))) {
					for (CSVRecord record : parser) {
						key = record.get(0);
						List<StreamingResource> resources = getResources(record.get(1));
						if (resources != null) {
							StreamingResourceListStream ll = (StreamingResourceListStream) resMap.get(key);
							if (ll == null) {
								ll = new StreamingResourceListStream();
								resMap.put(key, ll);
							}
							ll.getList().addAll(resources);
						} else if (getLogger().isWarnEnabled()) {
							getLogger().warn("Missing fastq or fasta resource: " + record.get(1));
						}
					}
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		}
		if (otherResources != null) {
			if (orKey == null) {
				getLogger().warn("Extra streaming resources key must not be null.");
			} else {
				if (resMap.containsKey(orKey)) {
					getLogger().warn("Extra streaming resource not under unique key: " + orKey);
				} else {
					resMap.put(orKey, otherResources);
				}
			}
		}
		return resMap;
	}

	protected List<StreamingResource> getResources(String pathOrURL) {
		try {
			return Collections.singletonList(new StreamingURLResource(new URL(pathOrURL)));
		} catch (MalformedURLException e) {
			List<File> fastqs = getProject().fastqFilesFromPath(pathOrURL, fastqType);
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
