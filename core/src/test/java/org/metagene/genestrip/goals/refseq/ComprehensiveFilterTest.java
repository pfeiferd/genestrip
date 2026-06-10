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
package org.metagene.genestrip.goals.refseq;

import org.apache.commons.io.FileUtils;
import org.junit.BeforeClass;
import org.junit.Test;
import org.metagene.genestrip.*;
import org.metagene.genestrip.goals.*;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.FastqKMerMatcher;
import org.metagene.genestrip.store.KMerSortedArray;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.net.URL;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertTrue;

public class ComprehensiveFilterTest extends ComprehensiveMatchTest {
    @BeforeClass()
    public static void clearDB() throws IOException {
        GSProject project = createProject("human_virus", null);
        GSMaker maker = new GSMaker(project);
        maker.getGoal(GSGoalKey.CLEAR).make();
        maker.dumpAll();
    }

    @Override
    public void testUpdate() throws IOException {
        // Just to avoid running the test from the superclass ...
    }

    @Override
    public void testKrakenOutput() throws IOException {
        // Just to avoid running the test from the superclass ...
    }

    // This test relies on the correctness of the matcher
    @Test
    public <P extends GSProject> void testFilterOutput() throws IOException {
        GSProject project = createProject("human_virus", null);
        project.initConfigParam(GSConfigKey.THREADS, 0);
        project.initConfigParam(GSConfigKey.WRITE_FILTERED_FASTQ, true);
        project.initConfigParam(GSConfigKey.GZIP_FASTQ_OUTPUT, false);
        project.logParamMap();

        createProjectGoal(project).make();

        GSMaker<P>  maker = new GSMaker<P>((P) project) {
            // A whole lot of overriding for testing purposes...
            protected MatchResultGoal<P> createGoalChainForMatchResult(boolean lr, String key, String... pathsOrURLs) {
                ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapGoal = new FastqMapGoal(getProject(), true,
                        getGoal(GSGoalKey.SETUP)) {
                    @Override
                    protected void doMakeThis() {
                        Map<String, StreamingResourceStream> map = createFastqMap(key, pathsOrURLs, null, null, null);
                        set(map);
                        if (getLogger().isInfoEnabled()) {
                            getLogger().info("Derived fastq map: " + map);
                        }
                    }
                };

                ObjectGoal<Map<String, StreamingResourceStream>, P> fastqMapTransfGoal = new FastqMapTransformGoal(
                        getProject(), true, fastqMapGoal, getGoal(GSGoalKey.SETUP));

                FastqDownloadsGoal<P> fastqDownloadsGoal = new FastqDownloadsGoal(getProject(), true, fastqMapGoal, fastqMapTransfGoal,
                        getGoal(GSGoalKey.SETUP));

                LoadDBGoal<P> loadDBGoal = getLoadDBGoal();

                return new MatchResultGoal<P>(getProject(), (lr ? GSGoalKey.MATCHRESLR : GSGoalKey.MATCHRES), fastqMapTransfGoal, loadDBGoal,
                        getExecutionContext(getProject()), getGoal(GSGoalKey.SETUP), fastqDownloadsGoal) {
                    @Override
                    protected FastqKMerMatcher createMatcher(KMerSortedArray<SmallTaxTree.SmallTaxIdNode> store, SmallTaxTree taxTree, ExecutionContext bundle, boolean withProbs, String dbMD5) {
                        return new FastqKMerMatcher(store, intConfigValue(GSConfigKey.INITIAL_READ_SIZE_BYTES),
                                intConfigValue(GSConfigKey.THREAD_QUEUE_SIZE), bundle, withProbs, intConfigValue(GSConfigKey.MAX_KMER_RES_COUNTS),
                                taxTree, intConfigValue(GSConfigKey.MAX_CLASSIFICATION_PATHS),
                                doubleConfigValue(GSConfigKey.MAX_READ_TAX_ERROR_COUNT),
                                doubleConfigValue(GSConfigKey.MAX_READ_CLASS_ERROR_COUNT),
                                booleanConfigValue(GSConfigKey.WRITE_ALL),
                                intConfigValue(GSConfigKey.MIN_KMERS_FOR_CLASS),
                                dbMD5) {
                            @Override
                            protected boolean isProgressBar() {
                                return booleanConfigValue(GSConfigKey.PROGRESS_BAR);
                            }

                            @Override
                            protected String getProgressBarTaskName() {
                                return getKey().getName();
                            }

                            // If even a single k-mer in the read is from a requested node,
                            // we must print that read for the comparison.
                            @Override
                            protected void rewriteInput(ReadEntry readStruct, OutputStream out) throws IOException {
                                MatcherReadEntry matcherReadEntry = (MatcherReadEntry) readStruct;
                                for (int i = 0; i < matcherReadEntry.usedPaths; i++) {
                                    if (matcherReadEntry.readTaxIdNode[i].isRequested()) {
                                        super.rewriteInput(readStruct, out);
                                        break;
                                    }
                                }
                            }
                        };
                    }
                };
            }

        };
        maker.getGoal(GSGoalKey.DB).make();
        maker.getGoal(GSGoalKey.INDEX).make();
        maker.dumpAll();

        File toFilter = new File(project.getFastaDir(), "test.fasta.gz");
        maker.match(false, false, "test", toFilter.toString());
        maker.filter(true, "test2", toFilter.toString());

        File file1 = new File(project.getFastqDir(), getProjectName() + "_matchres_test.fastq");
        File file2 = new File(project.getFastqDir(), getProjectName() + "_filter_test2.fastq");
        System.out.println("file1: " + file1);
        System.out.println("file2: " + file2);
        assertTrue(FileUtils.contentEquals(file1, file2));
        maker.dumpAll();
    }

    protected FileListGoal<GSProject> createProjectGoal(GSProject project) {
        return new HumanVirusProjectGoal(project);
    }

    protected String getProjectName() {
        return "human_virus";
    }

    public class HumanVirusProjectGoal extends FileListGoal<GSProject> {
        @SafeVarargs
        public HumanVirusProjectGoal(GSProject project, Goal<GSProject>... dependencies) {
            super(project, new GoalKey.DefaultGoalKey(getProjectName()), (List<File>) null, dependencies);
            addFile(new File(project.getProjectDir(), "taxids.txt"));
            addFile(new File(project.getProjectDir(), "categories.txt"));
            addFile(new File(project.getProjectDir(), "fasta/test.fasta.gz"));
        }

        @Override
        protected void makeFile(File file) throws IOException {
            if (!getProject().getCommon().getBaseDir().exists()) {
                getProject().getCommon().getBaseDir().mkdir();
            }
            if (!getProject().getProjectsDir().exists()) {
                getProject().getProjectsDir().mkdir();
            }
            if (!getProject().getProjectDir().exists()) {
                getProject().getProjectDir().mkdir();
            }
            if (!getProject().getFastaDir().exists()) {
                getProject().getFastaDir().mkdir();
            }
            if (!getProject().getFastqDir().exists()) {
                getProject().getFastqDir().mkdir();
            }
            if (!file.exists()) {
                ClassLoader classLoader = getClass().getClassLoader();
                URL url = classLoader.getResource("projects/" + getResFolderStr(file.getName()) + "/" + file.getName());
                if (url == null) {
                    url = classLoader.getResource("projects/viral/" + file.getName());
                }
                File orgFile = new File(url.getFile());
                if (!file.exists()) {
                    Files.copy(orgFile.toPath(), file.toPath());
                }
            }
        }

        protected String getResFolderStr(String fileName) {
            return getProjectName();
        }
    }
}
