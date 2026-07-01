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

/**
 * Renders a database's taxonomy tree as an indented, connected line drawing onto a
 * {@link Graphics2D}, labelling each node with its name, tax id and stored k-mer count.
 * Node indentation and dashing reflect the evolutionary distances computed by
 * {@link EvoDistanceEstimator}. Concrete subclasses supply the layout and rendering
 * parameters via {@link #initParams()}.
 */
public abstract class TaxTreePainter {
    private static final DecimalFormat DF = new DecimalFormat("0.000", new DecimalFormatSymbols(Locale.US));

    private final Database database;

    /** The font used to render standard (non-requested) node labels. */
    protected Font font;
    /** The font used to render labels of requested nodes. */
    protected Font requestedFont;
    /** The vertical spacing between lines, as a multiple of the font height. */
    protected double spacingFactor;
    /** The standard horizontal indent per level, as a multiple of the line height. */
    protected double indentFactor;
    /** The horizontal gap in pixels between a node's connector line and its text. */
    protected int xTextGap;
    /** The scaling factor applied to per-node distance- or k-mer-based indentation. */
    protected double nodeIndentFactor;
    /** Whether per-node indentation is derived from evolutionary distance rather than k-mer count. */
    protected boolean distanceIndent;
    /** The distance above which a node's connector line is drawn dashed. */
    protected double largeDistanceThreshold;
    /** Whether the longest evolutionary path is highlighted (in red). */
    protected boolean markLongestPath;
    /** Whether a node's taxonomic rank is included in its label. */
    protected boolean showRank;
    /** Whether a node's evolutionary distance is included in its label. */
    protected boolean showDistance;
    /** Whether a node's distance portion is included in its label. */
    protected boolean showDistancePortion;

    private BasicStroke dashed;

    /**
     * Creates a painter for the taxonomy tree of the given database.
     *
     * @param database the database whose taxonomy tree and stats are rendered
     */
    public TaxTreePainter(Database database) {
        this.database = database;
        initParams();
    }

    /**
     * Initializes the layout, font and rendering parameters. Called from the constructor.
     */
    protected abstract void initParams();

    /**
     * Paints the entire taxonomy tree onto the given graphics context.
     *
     * @param graphics the graphics context to draw onto
     * @return the size (width and height) required by the drawing
     */
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

    /**
     * Returns a progress bar to track painting, or {@code null} (the default) for none.
     *
     * @param update the progress update supplying current and maximum counts
     * @return a progress bar to track painting, or {@code null} for none
     */
    protected ProgressBar getProgressBar(GSProgressUpdate update) {
        return null;
    }

    /**
     * Builds the glyph vector for a node's label text in the appropriate font.
     *
     * @param graphics the graphics context supplying the font render context
     * @param distanceInfo the distance info identifying the node to label
     * @param stats the per-tax-id k-mer counts used in the label
     * @return the glyph vector for the node's label text
     */
    protected GlyphVector createGlyphVector(Graphics2D graphics, EvoDistanceEstimator.DistanceInfo distanceInfo, Object2LongMap<String> stats) {
        String text = getNodeText(distanceInfo, stats);
        return getFont(graphics, distanceInfo.getNode(), stats).createGlyphVector(graphics.getFontRenderContext(), text);
    }

    /**
     * Returns the font to use for the given node's label, using the requested font for
     * requested nodes and the standard font otherwise.
     *
     * @param graphics the graphics context
     * @param node the node whose label font is requested
     * @param stats the per-tax-id k-mer counts
     * @return the font to use for the node's label
     */
    protected Font getFont(Graphics2D graphics, SmallTaxTree.SmallTaxIdNode node, Object2LongMap<String> stats) {
        return node.isRequested() ? requestedFont : font;
    }

    /**
     * Returns the maximum node depth (level) occurring in the tree.
     *
     * @param tree the taxonomy tree to inspect
     * @return the maximum node level occurring in the tree
     */
    protected int getMaxLevel(SmallTaxTree tree) {
        int maxLevel = 0;
        Iterator<SmallTaxTree.SmallTaxIdNode> iterator = tree.iterator();
        while (iterator.hasNext()) {
            SmallTaxTree.SmallTaxIdNode node = iterator.next();
            maxLevel = Math.max(maxLevel, node.getLevel());
        }

        return maxLevel;
    }

    /**
     * Whether the given node is the last child of its parent (used to stop drawing the
     * parent's vertical connector line).
     *
     * @param node the node to test
     * @return {@code true} if the node is the last child of its parent (or has no parent)
     */
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

    /**
     * Returns the extra horizontal indent for a node, derived either from its distance
     * portion or from its k-mer count relative to the maximum, depending on configuration.
     *
     * @param info the distance info for the node
     * @param kmers the node's stored k-mer count
     * @param maxKMers the maximum k-mer count across all nodes
     * @return the extra horizontal indent in pixels for the node
     */
    protected int indentForNode(EvoDistanceEstimator.DistanceInfo info, long kmers, long maxKMers) {
        if (distanceIndent) {
            // distance == 1 gives no useful indentation (e.g. from nodes near the root)
            return info.getDistance() == 1 ? 0 : (int) (info.getDistancePortion() * nodeIndentFactor);
        } else {
            return (int) ((((double) kmers) / maxKMers) * nodeIndentFactor);
        }
    }

    /**
     * Returns the total horizontal offset of the given level, summing the per-node
     * indents up to that level plus the standard indent for each level.
     *
     * @param indentNodeWidths the per-level extra indent widths
     * @param level the level whose total offset is computed
     * @param stdIndentWidth the standard indent width per level
     * @return the total horizontal offset in pixels of the given level
     */
    protected int getTotalIndentWidth(int[] indentNodeWidths, int level, int stdIndentWidth) {
        int sum = 0;
        for (int i = 0; i <= level; i++) {
            sum += indentNodeWidths[i];
        }
        return sum + level * stdIndentWidth;
    }

    /**
     * Builds a node's label, comprising its name, tax id (optionally rank), stored k-mer
     * count and, if enabled, its distance and distance portion.
     *
     * @param distanceInfo the distance info identifying the node to label
     * @param stats the per-tax-id k-mer counts used in the label
     * @return the label text for the node
     */
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
