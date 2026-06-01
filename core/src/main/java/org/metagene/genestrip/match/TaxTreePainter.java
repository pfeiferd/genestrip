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
package org.metagene.genestrip.match;

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import me.tongfei.progressbar.ProgressBar;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.util.progressbar.GSProgressUpdate;

import java.awt.*;
import java.awt.font.GlyphVector;
import java.awt.geom.Rectangle2D;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Iterator;
import java.util.Locale;
import java.util.Map;
import java.util.function.BiConsumer;

public abstract class TaxTreePainter {
    private static final DecimalFormat DF = new DecimalFormat("0.000", new DecimalFormatSymbols(Locale.US));

    private final Database database;

    protected Font font;
    protected Font requestedFont;
    protected double spacingFactor;
    protected double indentFactor;
    protected int xTextGap;
    protected double nodeIndentFactor;
    protected boolean distanceIndent;
    protected double largeDistanceThreshold;
    protected boolean markLongestPath;
    protected boolean showRank;
    protected boolean showDistance;
    protected boolean showDistancePortion;

    private BasicStroke dashed;

    public TaxTreePainter(Database database) {
        this.database = database;
        initParams();
    }

    protected abstract void initParams();

    public Dimension paintTree(Graphics2D graphics) {
        graphics.setFont(font);
        int ascent = graphics.getFontMetrics().getAscent();
        int descent = graphics.getFontMetrics().getDescent();
        int fontHeight = ascent + descent;
        int lineHeight = (int) Math.round(spacingFactor * fontHeight);
        int stdIndentWidth = (int) Math.round(indentFactor * lineHeight);

        dashed = new BasicStroke(fontHeight / 10f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER,
                5.0f, new float[]{ 5.0f }, 0.0f);

        graphics.setColor(Color.black);
        graphics.setStroke(new BasicStroke(fontHeight / 10f));

        SmallTaxTree tree = database.getTaxTree();
        Object2LongMap<String> stats = database.getStats();
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

        int[] c = new int[1];
        Map<SmallTaxTree.SmallTaxIdNode, EvoDistanceEstimator.DistanceInfo> distances = new EvoDistanceEstimator().computeDistances(database);

        GSProgressUpdate update = new GSProgressUpdate() {
            @Override
            public long current() {
                return c[0];
            }

            @Override
            public long max() {
                // '- 1' to omit the null entry...
                return stats.size() - 1;
            }
        };
        try (ProgressBar pb = getProgressBar(update)) {
            Iterator<SmallTaxTree.SmallTaxIdNode> iterator = tree.iterator();
            for (c[0] = 0; iterator.hasNext(); c[0]++) {
                SmallTaxTree.SmallTaxIdNode node = iterator.next();
                EvoDistanceEstimator.DistanceInfo distanceInfo = distances.get(node);
                long s = stats.getOrDefault(node.getTaxId(), 0L);

                int level = node.getLevel();
                indentNodeWidths[level] = indentForNode(distanceInfo, s, maxKMers[0]);
                drawLines[level] = true;

                SmallTaxTree.SmallTaxIdNode parentBranch = null;
                if (markLongestPath) {
                    SmallTaxTree.SmallTaxIdNode parent = node.getParent();
                    if (parent != null) {
                        parentBranch = distances.get(parent).getBranch();
                    }
                }
                int xLineOffSet = (int) (lineHeight * 0.25);
                boolean lastChild = isLastChild(node);
                for (int j = 1; j <= level; j++) {
                    int xPos = getTotalIndentWidth(indentNodeWidths, j - 1, stdIndentWidth) + xLineOffSet;
                    int yTop = c[0] * lineHeight;
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
                        Color oldColor = graphics.getColor();
                        if (node == parentBranch) {
                            graphics.setColor(Color.RED);
                        }
                        if (distanceInfo.getDistance() > largeDistanceThreshold) {
                            Stroke oldStroke = graphics.getStroke();
                            graphics.setStroke(dashed);
                            graphics.drawLine(xPos, yBottom, xEnd, yBottom);
                            graphics.setStroke(oldStroke);
                        } else {
                            graphics.drawLine(xPos, yBottom, xEnd, yBottom);
                        }
                        graphics.setColor(oldColor);
                    }
                }
                if (lastChild) {
                    drawLines[level] = false;
                }

                // Paint the text
                GlyphVector gv = createGlyphVector(graphics, distanceInfo, stats);
                int x = getTotalIndentWidth(indentNodeWidths, level, stdIndentWidth);
                int y = c[0] * lineHeight + (lineHeight - fontHeight) / 2 + ascent;
                // Clear background behind text so that potential tree lines are out of the way.
                graphics.setColor(Color.WHITE);
                Rectangle2D r = gv.getVisualBounds(); // The visual bounds are not entirely correct - hence corrections below
                graphics.fillRect(x, y - (int) r.getHeight(), (int) r.getWidth() + fontHeight, (int) r.getHeight() + descent);
                graphics.setColor(Color.BLACK);
                graphics.drawGlyphVector(gv, x, y);
                maxWidth = Math.max(maxWidth, x + (int) r.getWidth());
            }

            return new Dimension(maxWidth + fontHeight, c[0] * lineHeight);
        }
    }

    protected ProgressBar getProgressBar(GSProgressUpdate update) {
        return null;
    }

    protected GlyphVector createGlyphVector(Graphics2D graphics, EvoDistanceEstimator.DistanceInfo distanceInfo, Object2LongMap<String> stats) {
        String text = getNodeText(distanceInfo, stats);
        return getFont(graphics, distanceInfo.getNode(), stats).createGlyphVector(graphics.getFontRenderContext(), text);
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

    protected int indentForNode(EvoDistanceEstimator.DistanceInfo info, long kmers, long maxKMers) {
        if (distanceIndent) {
            // distance == 1 gives no useful indentation (e.g. from nodes near the root)
            return info.getDistance() == 1 ? 0 : (int) (info.getDistancePortion() * nodeIndentFactor);
        } else {
            return (int) ((((double) kmers) / maxKMers) * nodeIndentFactor);
        }
    }

    protected int getTotalIndentWidth(int[] indentNodeWidths, int level, int stdIndentWidth) {
        int sum = 0;
        for (int i = 0; i <= level; i++) {
            sum += indentNodeWidths[i];
        }
        return sum + level * stdIndentWidth;
    }

    protected String getNodeText(EvoDistanceEstimator.DistanceInfo distanceInfo, Object2LongMap<String> stats) {
        StringBuilder sb = new StringBuilder();
        SmallTaxTree.SmallTaxIdNode node = distanceInfo.getNode();
        sb.append(node.getName());
        sb.append(" (");
        sb.append(node.getTaxId());
        if (showRank) {
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
        if (showDistance) {
            sb.append(",d=");
            sb.append(DF.format(distanceInfo.getDistance()));
        }
        if (showDistancePortion) {
            sb.append(",dp=");
            sb.append(DF.format(distanceInfo.getDistancePortion()));
        }
        sb.append(']');
        return sb.toString();
    }
}
