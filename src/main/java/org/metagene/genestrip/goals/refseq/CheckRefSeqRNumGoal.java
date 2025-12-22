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

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;

import java.io.*;
import java.net.URL;
import java.nio.file.Files;

public class CheckRefSeqRNumGoal extends ObjectGoal<CheckRefSeqRNumGoal.Result, GSProject> {
    public enum Result { CURRENT, OUTDATED, UNKNOWN }

    private final FileGoal<GSProject> releaseNumberGoal;

    public CheckRefSeqRNumGoal(GSProject project, FileGoal<GSProject> releaseNumberGoal, Goal<GSProject>... deps) {
        super(project, GSGoalKey.CHECK_REFSEQ_RNUM, append(deps, releaseNumberGoal));
        this.releaseNumberGoal = releaseNumberGoal;
    }

    @Override
    protected void doMakeThis() {
        try {
            String url = buildHttpURL();
            try (InputStream inputStream = new URL(url).openStream();
                 ByteArrayOutputStream out = new ByteArrayOutputStream()) {
                for (int c = inputStream.read(); c != -1; c = inputStream.read()) {
                    out.write(c);
                }
                String currentReleaseNumber = out.toString().trim();
                byte[] encoded = Files.readAllBytes(releaseNumberGoal.getFile().toPath());
                String releaseNumber = new String(encoded).trim();

                boolean equal = currentReleaseNumber.equals(releaseNumber);
                if (!equal) {
                    if (getLogger().isWarnEnabled()) {
                        getLogger().warn("You are working with an outdated RefSeq release. Your release: " + releaseNumber + ". Current release: " + currentReleaseNumber);
                    }
                }
                set(equal ? Result.CURRENT : Result.OUTDATED);
            }
        }
        catch (IOException e) {
            set(Result.UNKNOWN);
            if (getLogger().isErrorEnabled()) {
                getLogger().error("The current RefSeq release could not be determined.", e);
            }
        }
    }

    protected String buildHttpURL() {
        return getHttpBaseURL() + RefSeqDownloadGoal.RELEASE_FOLDER + "/" + RefSeqRNumDownloadGoal.RELEASE_NUMBER_FILE_NAME;
    }

    protected String getHttpBaseURL() {
        return stringConfigValue(GSConfigKey.REF_SEQ_HTTP_BASE_URL);
    }
}
