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
package org.metagene.genestrip.finertree.goals;

import it.unimi.dsi.fastutil.objects.Object2LongMap;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.finertree.FTConfigKey;
import org.metagene.genestrip.finertree.FTGoalKey;
import org.metagene.genestrip.finertree.FTProject;
import org.metagene.genestrip.finertree.cluster.DendrogramNode;
import org.metagene.genestrip.io.StreamProvider;
import org.metagene.genestrip.make.FileListGoal;
import org.metagene.genestrip.make.Goal;
import org.metagene.genestrip.make.ObjectGoal;
import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;

public class DengrogramLaTeXGoal<P extends FTProject> extends FileListGoal<P> {
    protected static final DecimalFormat LDF = new DecimalFormat("#,###", new DecimalFormatSymbols(Locale.US));
    private static final DecimalFormat DF = new DecimalFormat("0.000000", new DecimalFormatSymbols(Locale.US));
    private static final DecimalFormat DF2 = new DecimalFormat("0.#####", new DecimalFormatSymbols(Locale.US));

    private final ObjectGoal<Database, P> storeGoal;
    private final ObjectGoal<Map<SmallTaxTree.SmallTaxIdNode, DendrogramNode>, P> dendrogramGoal;
    private final Map<File, SmallTaxTree.SmallTaxIdNode> fileToNodeMap;
    private final double yScaleFactor ;
    private final double xScaleFactor;
    private final double tikzScale;
    private final boolean turn;
    private final boolean rescale;
    private Object2LongMap<String> stats;

    public DengrogramLaTeXGoal(P project, ObjectGoal<Database, P> storeGoal, ObjectGoal<Map<SmallTaxTree.SmallTaxIdNode, DendrogramNode>, P> dendrogramGoal, Goal<P>... deps) {
        super(project, FTGoalKey.DENDRO_LATEX, (List<File>) null, append(deps, storeGoal, dendrogramGoal));
        this.storeGoal = storeGoal;
        this.dendrogramGoal = dendrogramGoal;
        fileToNodeMap = new HashMap<>();

        yScaleFactor = doubleConfigValue(FTConfigKey.Y_FACTOR_LATEX);
        xScaleFactor = doubleConfigValue(FTConfigKey.X_FACTOR_LATEX);
        tikzScale = doubleConfigValue(FTConfigKey.TIKZ_SCALE_FACTOR);
        turn = booleanConfigValue(FTConfigKey.TURN_LATEX);
        rescale = booleanConfigValue(FTConfigKey.SIM_LOG_SCALING);
    }

    @Override
    public boolean isCleaned() {
        // Actual check would require make of other goals - not doable.
        return true;
    }

    @Override
    // Do not access kmerIntersectGoal here as it would trigger the related computation already...
    protected void provideFiles() {
        Database database = storeGoal.get();
        stats = database.getStats();
        SmallTaxTree tree = database.getTaxTree();
        Collection<SmallTaxTree.SmallTaxIdNode> parents = KMerIntersectCSVGoal.getNodesForPositions(database, (Collection<FTConfigKey.RefinementPosition>) configValue(FTConfigKey.REFINEMENT_POSITIONS));
        List<String> taxids = (List<String>) configValue(GSConfigKey.TAX_IDS);
        for (String taxid : taxids) {
            SmallTaxTree.SmallTaxIdNode parentNode = tree.getNodeByTaxId(taxid);
            if (parentNode != null) {
                parents.add(parentNode);
            }
        }
        for (SmallTaxTree.SmallTaxIdNode node : parents) {
            File matchFile = getProject().getOutputFile(getKey().getName(), node.getTaxId(), null, FTProject.FTFileType.TEX, false);
            addFile(matchFile);
            fileToNodeMap.put(matchFile, node);
        }
    }

    @Override
    protected void makeFile(File file) throws IOException {
        SmallTaxTree.SmallTaxIdNode parent = fileToNodeMap.get(file);
        DendrogramNode dendrogram = dendrogramGoal.get().get(parent);

        double[] minSim = new double[] { 1d };
        dendrogram.visit(new DendrogramNode.Visitor() {
            @Override
            public void preNode(DendrogramNode node) {
                if (node.getSimilarity() < minSim[0]) {
                    minSim[0] = node.getSimilarity();
                }
            }

            @Override
            public void postNode(DendrogramNode node) {
            }
        });
        double minLogSim = Math.log(minSim[0]);

        try (PrintStream out = new PrintStream(StreamProvider.getOutputStreamForFile(file))) {
            out.println("\\begin{figure}");
            out.println("\\begin{center}");
            out.print("\\begin{tikzpicture}[sloped,scale=");
            out.print(DF2.format(tikzScale));
            out.println("]");
            drawAxis(out, xScaleFactor, yScaleFactor, turn, minLogSim);
            drawDendrogram(out, parent, dendrogram, xScaleFactor, yScaleFactor, turn, minLogSim);
            out.println("\\end{tikzpicture}");
            out.print("\\caption{");
            out.print(parent.getName());
            out.print(" (");
            out.print(parent.getTaxId());
            out.print(")");
            if (parent.getRank() != null) {
                out.print(" [");
                out.print(parent.getRank().getName());
                out.print(" with ");
                out.print(LDF.format(stats.getOrDefault(parent.getTaxId(), 0)));
                out.print(" $k$-mers");
                out.print("]");
            }
            out.print("\\label{dendrogram");
            out.print(parent.getTaxId());
            out.println("}}");
            out.println("\\end{center}");
            out.println("\\end{figure}");
        }
    }

    protected void drawDendrogram(PrintStream out, SmallTaxTree.SmallTaxIdNode parent, DendrogramNode dendrogram, double xScaleFactor, double yScaleFactor, boolean turn, double minLogSim) {
        if (dendrogram == null) {
            return;
        }
        SmallTaxTree.SmallTaxIdNode[] children = parent.getSubNodes();
        int[] leafCounter = new int[1];
        int[] preCounter = new int[1];
        dendrogram.visit(new DendrogramNode.Visitor() {
            @Override
            public void preNode(DendrogramNode node) {
                // Exclude "OTHER" from display if it carries no information.
                if (node.getValueIndex() != children.length || node.getSimilarity() != 0) {
                    node.setValue(new IntDouble(preCounter[0], leafCounter[0]));
                }
                int index = node.getValueIndex();
                if (index >= 0) {
                    out.print("\\node ");
                    if (!turn) {
                        out.print("[rotate=90,anchor=east] ");
                    }
                    else {
                        out.print("[anchor=east] ");
                    }
                    out.print("(n");
                    out.print(preCounter[0]);
                    out.print(") at (");
                    if (turn) {
                        out.print("0,");
                    }
                    out.print(DF.format(xScaleFactor * leafCounter[0]));
                    if (!turn) {
                        out.print(",0");
                    }
                    out.print(") {");
                    if (index < children.length) {
                        out.print(children[index].getName());
                        out.print(" (");
                        out.print(children[index].getTaxId());
                        out.println(")");
                    }
                    else {
                        out.print("OTHER");
                    }
                    out.println("};");
                    leafCounter[0]++;
                }
                preCounter[0]++;
            }

            public void postNode(DendrogramNode node) {
                if (node.getValue() == null) {
                    return;
                }
                if (node.getValueIndex() == -1) {
                    double xPos = (((IntDouble) node.getChild1().getValue()).d + ((IntDouble) node.getChild2().getValue()).d) / 2;
                    IntDouble value = (IntDouble) node.getValue();
                    value.d = xPos;
                    out.print("\\node (n");
                    out.print(value.i);
                    out.print(") at (");
                    if (turn) {
                        out.print(DF.format(yScaleFactor * (1 - rescaleSim(node.getSimilarity(), minLogSim))));
                        out.print(",");
                        out.print(DF.format(xScaleFactor * xPos));
                    }
                    else {
                        out.print(DF.format(xScaleFactor * xPos));
                        out.print(",");
                        out.print(DF.format(yScaleFactor * (1 - rescaleSim(node.getSimilarity(), minLogSim))));
                    }
                    out.println(") {};");
                }
            }
        });
        preCounter[0] = 0;
        dendrogram.visit(new DendrogramNode.Visitor() {
            @Override
            public void preNode(DendrogramNode node) {
                if (node.getValueIndex() == -1) {
                    if (turn) {
                        drawEdge(node, node.getChild1(), turn);
                        drawEdge(node, node.getChild2(), turn);
                    }
                    else {
                        drawEdge(node.getChild1(), node, turn);
                        drawEdge(node.getChild2(), node, turn);
                    }
                }
            }

            protected void drawEdge(DendrogramNode from, DendrogramNode to, boolean turn) {
                if (from.getValue() == null || to.getValue() == null) {
                    return;
                }

                out.print("\\draw  (n");
                out.print(((IntDouble) from.getValue()).i);
                if (from.getValueIndex() == -1) {
                    out.print(".center");
                }
                else if (turn) {
                    out.print(".east");
                }
                out.print(") |- (n");
                out.print(((IntDouble) to.getValue()).i);
                if (to.getValueIndex() == -1) {
                    out.print(".center");
                }
                else if (turn) {
                    out.print(".east");
                }
                out.println(");");
            }

            @Override
            public void postNode(DendrogramNode node) {
            }
        });
    }

    protected double rescaleSim(double similarity, double minLogSim) {
        return rescale ? (1 - Math.log(similarity) / minLogSim) : similarity;
    }

    protected void drawAxis(PrintStream out, double xScaleFactor, double yScaleFactor, boolean turn, double minLogSim) {
        double xPos = -xScaleFactor -1.5;
        double yPos = yScaleFactor;
        out.print("\\draw[<-] (");
        out.print(turn ? "0" : DF.format(xPos));
        out.print(",");
        out.print(turn ? DF.format(xPos) : "0");
        out.print(") -- node[above]{Similarity} (");
        out.print(DF.format(turn ? yPos: xPos));
        out.print(",");
        out.print(DF.format(turn ? xPos : yPos));
        out.println(");");

        xPos = -xScaleFactor;
        out.print("\\draw (");
        out.print(turn ? "0" : DF.format(xPos));
        out.print(",");
        out.print(turn ? DF.format(xPos) : "0");
        out.print(") -- (");
        out.print(DF.format(turn ? yPos : xPos));
        out.print(",");
        out.print(DF.format(turn ? xPos : yPos));
        out.println(");");

        int max = 5;
        int expStep = (int) (- minLogSim / Math.log(10) / (max - 1));
        if (expStep <= 0) {
            expStep = 1;
        }
        for (int i = max; i >= (rescale ? 1 : 0); i--) {
            double xPosLeft = -xScaleFactor - 0.1;
            double v = rescale ? Math.pow(10, expStep * (i - max)) : ((double) i) / max;
            yPos = yScaleFactor * (1 - rescaleSim(v, minLogSim));
            if (yPos <= yScaleFactor) {
                out.print("\\draw (");
                out.print(DF.format(turn ? yPos : xPos));
                out.print(",");
                out.print(DF.format(turn ? xPos : yPos));
                out.print(") -- (");
                out.print(DF.format(turn ? yPos : xPosLeft));
                out.print(",");
                out.print(DF.format(turn ? xPosLeft : yPos));
                out.println(");");

                xPosLeft = -xScaleFactor - (turn ? 0.4 : 0.1);
                out.print(turn ? "\\node at (" : "\\node[left] at (");
                out.print(DF.format(turn ? yPos : xPosLeft));
                out.print(",");
                out.print(DF.format(turn ? xPosLeft : yPos));
                out.print(") {$");
                out.print(DF2.format(v));
                out.println("$};");
            }
        }
    }

    private static class IntDouble {
        public IntDouble(int i, double d) {
            this.i = i;
            this.d = d;
        }

        public int i;
        public double d;
    }
}