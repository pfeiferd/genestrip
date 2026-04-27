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
import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.GSProject;
import org.metagene.genestrip.make.FileGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import java.awt.*;
import java.awt.font.GlyphVector;
import java.awt.geom.Rectangle2D;
import java.io.*;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.function.BiConsumer;

public class SVGTaxTreeGoal<P extends GSProject> extends FileGoal<P> {
    private final Font font;
    private final Font requestedFont;
    private final double spacingFactor;
    private final double indentFactor;
    private final int xTextGap;
    private double nodeIndentFactor;
    private boolean logBasedIndent;
    private boolean addRank;
    private final ObjectGoal<Database, P> storeGoal;

    public SVGTaxTreeGoal(P project, GoalKey key, ObjectGoal<Database, P> storeGoal, Goal<P>... deps) {
        super(project, key, Goal.append(deps, storeGoal));
        this.storeGoal = storeGoal;
        int fontSize = intConfigValue(GSConfigKey.SVG_FONT_SIZE);
        font = new Font(stringConfigValue(GSConfigKey.SVG_FONT), Font.PLAIN, fontSize);
        requestedFont = booleanConfigValue(GSConfigKey.SVG_REQ_NODES) ? font.deriveFont(Font.BOLD) : font;
        spacingFactor = doubleConfigValue(GSConfigKey.SVG_LINE_HEIGHT_FACTOR);
        indentFactor = doubleConfigValue(GSConfigKey.SVG_INDENTATION_FACTOR);
        xTextGap = (int) (fontSize * doubleConfigValue(GSConfigKey.SVG_TEXT_GAP_FACTOR));
        nodeIndentFactor = doubleConfigValue(GSConfigKey.SVG_KMER_NODE_INDENT_FACTOR);
        logBasedIndent = booleanConfigValue(GSConfigKey.SVG_LOG_BASED_INDENT);
        addRank = booleanConfigValue(GSConfigKey.SVG_SHOW_RANK);
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
        Dimension d = paintTree(svgGenerator);
        svgGenerator.setSVGCanvasSize(d);

        // Finally, stream out SVG to the standard output using
        // UTF-8 encoding.
        try (Writer out = new OutputStreamWriter(new FileOutputStream(file), "UTF-8")) {
            svgGenerator.stream(out, false);
        }
    }

    protected Dimension paintTree(Graphics2D graphics) {
        graphics.setFont(font);
        int ascent = graphics.getFontMetrics().getAscent();
        int descent = graphics.getFontMetrics().getDescent();
        int fontHeight = ascent + descent;
        int lineHeight = (int) Math.round(spacingFactor * fontHeight);
        int stdIndentWidth = (int) Math.round(indentFactor * lineHeight);

        graphics.setColor(Color.black);
        graphics.setStroke(new BasicStroke(fontHeight / 10f));

        SmallTaxTree tree = storeGoal.get().getTaxTree();
        Iterator<SmallTaxTree.SmallTaxIdNode> iterator = tree.iterator();
        Object2LongMap<String> stats = storeGoal.get().getStats();
        long[] maxKMers = new long[1];
        stats.forEach(new BiConsumer<String, Long>() {
            @Override
            public void accept(String s, Long aLong) {
                if (s != null) {
                    maxKMers[0] = Math.max(maxKMers[0], aLong);
                }
            }
        });

        int maxWidth = 0;
        int maxLevel = getMaxLevel(tree);
        int[] indentNodeWidths = new int[maxLevel + 1];
        boolean[] drawLines = new boolean[indentNodeWidths.length];

        int c;
        int level = -1;
        for (c = 0; iterator.hasNext(); c++) {
            SmallTaxTree.SmallTaxIdNode node = iterator.next();

            level = node.getLevel();
            indentNodeWidths[level] = getIndentForNode(node, stats, maxKMers[0]);
            drawLines[level] = true;

            int xLineOffSet = (int) (lineHeight * 0.25);
            boolean lastChild = isLastChild(node);
            for (int j = 1; j <= level; j++) {
                int xPos = getTotalIndentWidth(indentNodeWidths, j - 1, stdIndentWidth) + xLineOffSet;
                int yTop = c * lineHeight;
                // Vertical lines
                if (drawLines[j]) {
                    int yBottom = yTop + lineHeight;
                    if (j == level && lastChild) {
                        yBottom -= lineHeight / 2;
                    }
                    graphics.drawLine(xPos, yTop, xPos, yBottom);
                }
                // Horizontal line
                if (j == level) {
                    int xEnd = getTotalIndentWidth(indentNodeWidths, j, stdIndentWidth) - xTextGap;
                    int yBottom = yTop + lineHeight / 2;
                    graphics.drawLine(xPos, yBottom, xEnd, yBottom);
                }
            }
            if (lastChild) {
                drawLines[level] = false;
            }

            // Paint the text
            GlyphVector gv = createGlyphVector(graphics, node, stats);
            int x = getTotalIndentWidth(indentNodeWidths, level, stdIndentWidth);
            int y = c * lineHeight + (lineHeight - fontHeight) / 2 + ascent;
            // Clear background behind text so that potential tree lines are out of the way.
            graphics.setColor(Color.WHITE);
            Rectangle2D r = gv.getVisualBounds(); // The visual bounds are not entirely correct - hence corrections below
            graphics.fillRect(x, y - (int) r.getHeight(), (int) r.getWidth() + fontHeight, (int) r.getHeight() + descent);
            graphics.setColor(Color.BLACK);
            graphics.drawGlyphVector(gv, x, y);
            maxWidth = Math.max(maxWidth, x + (int) r.getWidth());
        }
        return new Dimension(maxWidth + fontHeight, c * lineHeight);
    }

    protected GlyphVector createGlyphVector(Graphics2D graphics, SmallTaxTree.SmallTaxIdNode node, Object2LongMap<String> stats) {
        String text = getNodeText(node, stats);
        return getFont(graphics, node, stats).createGlyphVector(graphics.getFontRenderContext(), text);
    }

    protected Font getFont(Graphics2D graphics, SmallTaxTree.SmallTaxIdNode node, Object2LongMap<String> stats) {
        return node.isRequested() ? requestedFont : font;
    }

    protected int getMaxLevel(SmallTaxTree tree) {
        int maxLevel = 0;
        Iterator<SmallTaxTree.SmallTaxIdNode> iterator = tree.iterator();
        while (iterator.hasNext()) {
            SmallTaxTree.SmallTaxIdNode node = iterator.next();
            maxLevel = Math.max(maxLevel, node.getLevel());
        }

        return maxLevel;
    }

    protected boolean isLastChild(SmallTaxTree.SmallTaxIdNode node) {
        SmallTaxTree.SmallTaxIdNode parent = node.getParent();
        if (parent == null) {
            return true;
        }
        SmallTaxTree.SmallTaxIdNode[] children = parent.getSubNodes();
        for (int i = 0; i < children.length; i++) {
            if (children[i] == node && i == children.length - 1) {
                return true;
            }
        }
        return false;
    }

    protected int getIndentForNode(SmallTaxTree.SmallTaxIdNode node, Object2LongMap<String> stats, long maxKMers) {
        long s = stats.getOrDefault(node.getTaxId(), 0L);
        if (logBasedIndent) {
            double base = Math.log(maxKMers + 1);
            return (int) (Math.log(s + 1) / base * nodeIndentFactor);
        }
        else {
            return (int) ((((double) s) / maxKMers) * nodeIndentFactor);
        }
    }

    protected int getTotalIndentWidth(int[] indentNodeWidths, int level, int stdIndentWidth) {
        int sum = 0;
        for (int i = 0; i <= level; i++) {
            sum += indentNodeWidths[i];
        }
        return sum + level * stdIndentWidth;
    }

    protected String getNodeText(SmallTaxTree.SmallTaxIdNode node, Object2LongMap<String> stats) {
        StringBuilder sb = new StringBuilder();
        sb.append(node.getName());
        sb.append(" (");
        sb.append(node.getTaxId());
        if (addRank) {
            Rank r = node.getRank();
            if (r != null) {
                sb.append(", ");
                sb.append(r.getName());
            }
        }
        sb.append(')');
        sb.append(" [");
        long s = stats.getOrDefault(node.getTaxId(), 0L);
        sb.append(s);
        sb.append(']');
        return sb.toString();
    }
}
