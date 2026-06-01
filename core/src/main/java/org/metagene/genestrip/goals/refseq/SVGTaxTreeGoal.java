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

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import me.tongfei.progressbar.ProgressBar;
import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.match.TaxTreePainter;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.util.progressbar.GSProgressBarCreator;
import org.metagene.genestrip.util.progressbar.GSProgressUpdate;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import java.awt.*;
import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.List;

public class SVGTaxTreeGoal<P extends GSProject> extends FileGoal<P> {
    private final ObjectGoal<Database, P> storeGoal;

    public SVGTaxTreeGoal(P project, GoalKey key, ObjectGoal<Database, P> storeGoal, Goal<P>... deps) {
        super(project, key, Goal.append(deps, storeGoal));
        this.storeGoal = storeGoal;
    }

    @Override
    public List<File> getFiles() {
        return Collections.singletonList(getProject().getOutputFile(getKey().getName(), GSProject.GSFileType.SVG, false));
    }

    @Override
    protected void makeFile(File file) throws IOException {
        // Get a DOMImplementation.
        DOMImplementation domImpl =
                GenericDOMImplementation.getDOMImplementation();

        // Create an instance of org.w3c.dom.Document.
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);

        // Create an instance of the SVG Generator.
        SVGGraphics2D svgGenerator = new SVGGraphics2D(document);
        Dimension d = new TaxTreePainter(storeGoal.get()) {
            @Override
            protected void initParams() {
                int fontSize = intConfigValue(GSConfigKey.SVG_FONT_SIZE);
                font = new Font(stringConfigValue(GSConfigKey.SVG_FONT), Font.PLAIN, fontSize);
                requestedFont = booleanConfigValue(GSConfigKey.SVG_REQ_NODES) ? font.deriveFont(Font.BOLD) : font;
                spacingFactor = doubleConfigValue(GSConfigKey.SVG_LINE_HEIGHT_FACTOR);
                indentFactor = doubleConfigValue(GSConfigKey.SVG_INDENTATION_FACTOR);
                xTextGap = (int) (fontSize * doubleConfigValue(GSConfigKey.SVG_TEXT_GAP_FACTOR));
                nodeIndentFactor = doubleConfigValue(GSConfigKey.SVG_KMER_NODE_INDENT_FACTOR);
                distanceIndent = booleanConfigValue(GSConfigKey.SVG_DISTANCE_INDENT);
                markLongestPath = booleanConfigValue(GSConfigKey.SVG_MARK_LONGEST_PATH);
                showRank = booleanConfigValue(GSConfigKey.SVG_SHOW_RANK);
                showDistance = booleanConfigValue(GSConfigKey.SVG_SHOW_DISTANCE);
                showDistancePortion = booleanConfigValue(GSConfigKey.SVG_SHOW_DISTANCE_PORTION);
                largeDistanceThreshold = doubleConfigValue(GSConfigKey.SVG_TOO_LARGE_DISTANCE);
            }

            @Override
            protected ProgressBar getProgressBar(GSProgressUpdate update) {
                return booleanConfigValue(GSConfigKey.PROGRESS_BAR) ? GSProgressBarCreator.newGSProgressBar(getKey().getName(), " nodes", update, null) : null;
            }
        }.paintTree(svgGenerator);
        svgGenerator.setSVGCanvasSize(d);

        // Finally, stream out SVG to the standard output using
        // UTF-8 encoding.
        try (Writer out = new OutputStreamWriter(new FileOutputStream(file), "UTF-8")) {
            svgGenerator.stream(out, false);
        }
    }
}
