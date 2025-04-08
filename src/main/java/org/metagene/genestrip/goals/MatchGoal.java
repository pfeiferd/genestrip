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

import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.FileType;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.io.StreamingResourceStream;
import org.metagene.genestrip.make.*;
import org.metagene.genestrip.match.MatchingResult;
import org.metagene.genestrip.match.ResultReporter;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MatchGoal extends FileListGoal<GSProject> {
    private final ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapGoal;
    private final ObjectGoal<Map<String, MatchingResult>, GSProject> matchResGoal;
    private final Map<File, String> fileToKeyMap;
    private ResultReporter reporter;

    @SafeVarargs
    public MatchGoal(GSProject project, GoalKey key, ObjectGoal<Map<String, StreamingResourceStream>, GSProject> fastqMapGoal, ObjectGoal<Map<String, MatchingResult>, GSProject> matchResGoal, Goal<GSProject>... deps) {
        super(project, key, (List<File>) null, append(deps, fastqMapGoal, matchResGoal));
        this.fastqMapGoal = fastqMapGoal;
        this.matchResGoal = matchResGoal;
        fileToKeyMap = new HashMap<>();
    }

    @Override
    protected void provideFiles() {
        for (String key : fastqMapGoal.get().keySet()) {
            File matchFile = getProject().getOutputFile(getKey().getName(), key, null, FileType.CSV, false);
            addFile(matchFile);
            fileToKeyMap.put(matchFile, key);
        }
    }

    @Override
    protected void makeFile(File file) throws IOException {
        MatchingResult result = matchResGoal.get().get(fileToKeyMap.get(file));
        try (PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file))) {
            if (reporter == null) {
                reporter = new ResultReporter();
            }
            reporter.printMatchResult(result, out);
        }
    }
}
