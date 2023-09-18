package org.metagene.genestrip.goals.refseq;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileDownloadGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.util.ArraysUtil;
import org.metagene.genestrip.util.StreamProvider;

public class RefSeqFnaFilesDownloadGoal extends FileDownloadGoal<GSProject> {
	private static final CSVFormat FORMAT = CSVFormat.DEFAULT.builder().setDelimiter('\t').setRecordSeparator('\n')
			.build();

	private final RefSeqCatalogDownloadGoal catalogDLGoal;
	private final Collection<RefSeqCategory> categories;
	private List<File> files;
	private Map<File, RefSeqCategory> file2Cat;

	@SafeVarargs
	public RefSeqFnaFilesDownloadGoal(GSProject project, String name, Collection<RefSeqCategory> categories,
			RefSeqCatalogDownloadGoal catalogDLGoal, Goal<GSProject>... deps) {
		super(project, name, ArraysUtil.append(deps, catalogDLGoal));

		this.categories = Collections.unmodifiableCollection(categories);
		this.catalogDLGoal = catalogDLGoal;
	}

	public Collection<RefSeqCategory> getCategories() {
		return categories;
	}

	public RefSeqCategory getCategoryForFile(File file) {
		return file2Cat.get(file);
	}

	protected boolean isRelevantFileName(String filename) {
		return filename.endsWith(".genomic.fna.gz") || filename.endsWith(".genomic.fna");
	}

	protected RefSeqCategory getCategoryForFileName(String filename) {
		for (RefSeqCategory cat : categories) {
			if (filename.startsWith(cat.name() + ".")) {
				return cat;
			}
		}
		return null;
	}

	@Override
	protected String getFTPDir(File file) {
		return RefSeqCatalogDownloadGoal.RELEASE_FOLDER + "/" + file2Cat.get(file).getDirectory();
	}

	@Override
	public List<File> getFiles() {
		if (files == null) {
			files = new ArrayList<File>();
			file2Cat = new HashMap<File, RefSeqCategory>();

			File installedFiles = catalogDLGoal.getInstalledFilesFile();

			try {
				Reader in = new InputStreamReader(StreamProvider.getInputStreamForFile(installedFiles));
				Iterable<CSVRecord> records = FORMAT.parse(in);
				for (CSVRecord record : records) {
					String filename = record.get(1);
					RefSeqCategory cat = getCategoryForFileName(filename);
					if (cat != null) {
						if (isRelevantFileName(filename)) {
							File file = new File(getProject().getConfig().getRefSeqDir(), filename);
							files.add(file);
							file2Cat.put(file, cat);
						}
					}
				}
				in.close();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
		return files;
	}
}
