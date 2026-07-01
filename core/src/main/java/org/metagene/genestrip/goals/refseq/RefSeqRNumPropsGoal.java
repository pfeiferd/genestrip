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

import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Properties;

/**
 * Goal that writes a {@code release.properties} file recording the RefSeq release number the
 * database was built from.
 *
 * @param <P> the project type
 */
public class RefSeqRNumPropsGoal<P extends GSProject> extends FileListGoal<P> {
    /** Name of the properties file recording the RefSeq release number. */
    public static final String RELEASE_PROPS_FILE_NAME = "release.properties";

    private final ObjectGoal<CheckRefSeqRNumGoal.Result, P> checkRefSeqRNumGoal;

    /**
     * Creates the goal, wiring the release-check goal whose side effect records the release number.
     *
     * @param project   the project this goal belongs to
     * @param checkRNum the release-check goal whose side effect records the release number
     * @param deps      additional goals this goal depends on
     */
    @SafeVarargs
    public RefSeqRNumPropsGoal(P project, ObjectGoal<CheckRefSeqRNumGoal.Result, P> checkRNum, Goal<P>... deps) {
        super(project, GSGoalKey.REFSEQ_PROP, new File(project.getCommon().getRefSeqDir(), RELEASE_PROPS_FILE_NAME), deps);
        this.checkRefSeqRNumGoal = checkRNum;
    }

    @Override
    protected void makeFile(File storeFile) throws IOException {
        // Ensure side effect of pickung up the release number in the props happens:
        checkRefSeqRNumGoal.make();
        String release = getProject().getAdditionalProperty(GSProject.REFSEQ_RELEASE);
        Properties props = new Properties();
        if (release != null) {
            props.setProperty(GSProject.REFSEQ_RELEASE, release);
        }
        try (FileWriter writer = new FileWriter(storeFile)) {
            props.store(writer, null);
        }
    }
}