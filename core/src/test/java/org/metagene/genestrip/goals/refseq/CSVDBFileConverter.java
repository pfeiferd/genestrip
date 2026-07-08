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

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.metagene.genestrip.GSCommon;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSGoalKey;
import org.metagene.genestrip.GSMaker;
import org.metagene.genestrip.GSProject;

import java.io.*;

public class CSVDBFileConverter {
    protected static final CSVFormat CSV_FORMAT = CSVFormat.DEFAULT.builder().setQuote(null).setCommentMarker('#')
            .setDelimiter(';').setRecordSeparator('\n').build();

    private final File baseDir;

    public CSVDBFileConverter(File file) {
        this.baseDir = file;
    }

    // This format also works for KrakenUniq to build the custom database
    public void csv2KuMapFileFormat(String db) throws IOException {
        GSCommon config = new GSCommon(baseDir);
        GSProject project = new GSProject(config, db, null, null, null, null, null, null,
                null, null, null, false);
        project.initConfigParam(GSConfigKey.THREADS, -1);

        GSMaker maker = new GSMaker(project);

        ExtractRefSeqCSVGoal extractRefSeqCSVGoal = (ExtractRefSeqCSVGoal) maker.getGoal(GSGoalKey.EXTRACT_REFSEQ_CSV);
        //extractRefSeqCSVGoal.make();

        try (CSVParser parser = CSV_FORMAT
                .parse(new InputStreamReader(new FileInputStream(extractRefSeqCSVGoal.getFile())))) {
            File kuInputFile = new File(project.getResultsDir(), db + "_ku.map");
            try (PrintStream out = new PrintStream(new FileOutputStream(kuInputFile))) {
                int i = 0;
                for (CSVRecord record : parser.getRecords()) {
                    if (i > 0) {
                        String descr = record.get(0);
                        String taxid = record.get(1);
                        out.print(descr);
                        out.print("|kraken:taxid|");
                        out.print(taxid);
                        out.print('\t');
                        out.print(taxid);
                        out.println();
                    }
                    i++;
                }
            }
        }
    }

    public static void main(String[] args) throws IOException {
        CSVDBFileConverter converter = new CSVDBFileConverter(new File("./data"));
        converter.csv2KuMapFileFormat(args[0]);
    }
}
