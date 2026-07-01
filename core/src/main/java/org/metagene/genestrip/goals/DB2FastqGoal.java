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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.GSProject.GSFileType;
import org.metagene.genestrip.fastqgen.KMerFastqGenerator;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

/**
 * Generates FASTQ files from the k-mers stored in the database, writing one output file per configured tax id
 * (or a single {@code total} file when none are given). A trailing {@code +} on a tax id enables read descriptions.
 *
 * @param <P> the project type
 */
public class DB2FastqGoal<P extends GSProject> extends FileListGoal<P> {
    /** The output name used when no specific tax ids are configured. */
    public static final String TOTAL = "total";

    private final ObjectGoal<Set<SmallTaxIdNode>, P> taxNodesGoal;
    private final ObjectGoal<Database, P> storeGoal;
    private final Map<File, String> fileToTaxid;

    private final String[] taxids;
    private final boolean[] withDescs;

    /**
     * Creates the goal for generating FASTQ files from the database's stored k-mers.
     *
     * @param project the project this goal belongs to
     * @param key the goal key identifying this goal
     * @param taxNodesGoal the goal supplying the set of taxonomy nodes
     * @param storeGoal the goal supplying the database to read k-mers from
     * @param deps additional goals this goal depends on
     */
    @SafeVarargs
    public DB2FastqGoal(P project, GoalKey key, ObjectGoal<Set<SmallTaxIdNode>, P> taxNodesGoal,
                        ObjectGoal<Database, P> storeGoal, Goal<P>... deps) {
        super(project, key, (List<File>) null, true, Goal.append(deps, taxNodesGoal, storeGoal));
        this.taxNodesGoal = taxNodesGoal;
        this.storeGoal = storeGoal;
        this.fileToTaxid = new HashMap<>();

        List<String> taxidList = (List<String>) getProject().configValue(GSConfigKey.TAX_IDS);
        if (taxidList != null && !taxidList.isEmpty()) {
            taxids = taxidList.toArray(new String[taxidList.size()]);
            withDescs = new boolean[taxids.length];
            for (int i = 0; i < taxids.length; i++) {
                if (taxids[i].endsWith("+")) {
                    taxids[i] = taxids[i].substring(0, taxids[i].length() - 1);
                    withDescs[i] = true;
                }
            }
        } else {
            taxids = new String[]{ TOTAL };
            withDescs = new boolean[]{true};
        }
    }

    @Override
    protected void provideFiles() {
        for (String taxid : taxids) {
            File file = getOutputFile(taxid);
            addFile(file);
            fileToTaxid.put(file, taxid);
        }
    }

    /**
     * Returns the output FASTQ file for the given tax id.
     *
     * @param taxid the tax id (or {@link #TOTAL}) the file is generated for
     * @return the output FASTQ file for the given tax id
     */
    protected File getOutputFile(String taxid) {
        return getProject().getOutputFile(getKey().getName(), taxid, null, GSFileType.FASTQ_RES, booleanConfigValue(GSConfigKey.GZIP_FASTQ_OUTPUT));
    }

    @Override
    protected void makeFile(File file) throws IOException {
        KMerFastqGenerator generator = new KMerFastqGenerator(storeGoal.get());
        String taxid = fileToTaxid.get(file);
        boolean withDesc = false;
        for (int i = 0; i < taxids.length; i++) {
            if (taxids[i].equals(taxid)) {
                withDesc = withDescs[i];
                break;
            }
        }
        generator.generateFastq(file, TOTAL.equals(taxid) ? null : taxid, getProject().getName() + ":", withDesc);
    }
}
