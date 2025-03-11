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

import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;

import org.metagene.genestrip.store.Database;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.SmallTaxTree;
import org.metagene.genestrip.tax.SmallTaxTree.SmallTaxIdNode;

import it.unimi.dsi.fastutil.objects.Object2LongMap;

public class ResultReporter {
    private static final DecimalFormat DF = new DecimalFormat("0.00000000", new DecimalFormatSymbols(Locale.US));

    public ResultReporter() {
    }

    public void printStoreInfo(Database database, PrintStream out) {
        Object2LongMap<String> stats = database.getStats();

        out.println("name;rank;taxid;stored kmers;");

        out.print("TOTAL;");
        out.print(Rank.NO_RANK);
        out.print(';');
        out.print("1;");
        out.print(stats.getLong(null));
        out.println(';');

        List<String> sortedTaxIds = new ArrayList<>(stats.keySet());
        SmallTaxTree taxTree = database.getTaxTree();
        taxTree.sortTaxidsViaTree(sortedTaxIds);

        for (String taxId : sortedTaxIds) {
            if (taxId != null) {
                SmallTaxIdNode taxNode = taxTree.getNodeByTaxId(taxId);
                if (taxNode != null) {
                    out.print(taxNode.getName());
                    out.print(';');
                    out.print(taxNode.getRank());
                    out.print(';');
                    out.print(taxNode.getTaxId());
                    out.print(';');
                    out.print(stats.getLong(taxId));
                    out.println(';');
                }
            }
        }
    }

    public static void printMDColumnInfo(PrintStream ps) {
        ps.print('|');
        ps.print("Name");
        ps.print('|');
        ps.print("Description");
        ps.print('|');
        ps.println();

        ps.print('|');
        ps.print('-');
        ps.print('|');
        ps.print('-');
        ps.print('|');
        ps.println();

        for (MethodAndDescription methodAndDescription : getSortedMethodAndDescriptions()) {
            ps.print('|');
            ps.print('`');
            ps.print(methodAndDescription.getDescription().name());
            ps.print('`');
            ps.print('|');
            ps.print(methodAndDescription.getDescription().value());
            ps.print('|');
            ps.println();
        }
    }

    private static List<MethodAndDescription> getSortedMethodAndDescriptions() {
        List<MethodAndDescription> methodAndDescriptions = new ArrayList<>();
        Method[] methods = CountsPerTaxid.class.getMethods();
        for (Method method : methods) {
            MDCDescription description = method.getAnnotation(MDCDescription.class);
            if (description != null) {
                methodAndDescriptions.add(new MethodAndDescription(method, description));
            }
        }
        Collections.sort(methodAndDescriptions);
        return methodAndDescriptions;
    }

    public void printMatchResult(MatchingResult res, PrintStream out) {
        List<MethodAndDescription> methodAndDescriptions = getSortedMethodAndDescriptions();

        for (MethodAndDescription methodAndDescription : methodAndDescriptions) {
            if (methodAndDescription.getDescription().pos() != 1001 || res.isWithMaxKMerCounts()) {
                out.print(methodAndDescription.getDescription().name());
                out.print(';');
            }
        }
        out.println();

        List<CountsPerTaxid> values = new ArrayList<>(res.getTaxid2Stats().values());
        Collections.sort(values);

        for (CountsPerTaxid counts : values) {
            for (MethodAndDescription methodAndDescription : methodAndDescriptions) {
                if (methodAndDescription.getDescription().pos() != 1001 || res.isWithMaxKMerCounts()) {
                    Method method = methodAndDescription.getMethod();
                    try {
                        Object value = method.invoke(counts);
                        if (value instanceof Double) {
                            double v = ((Double) value).doubleValue();
                            if (!Double.isNaN(v) && !Double.isInfinite(v)) {
                                out.print(v);
                            }
                        } else if (methodAndDescription.getDescription().pos() == 1001 && res.isWithMaxKMerCounts()) {
                            short[] maxKMerCounts = counts.getMaxKMerCounts();
                            if (maxKMerCounts != null) {
                                for (int i = 0; i < maxKMerCounts.length; i++) {
                                    if (i > 0) {
                                        out.print(';');
                                    }
                                    out.print(maxKMerCounts[i]);
                                }
                            }
                        } else {
                            out.print(value);
                        }
                        out.print(';');
                    } catch (IllegalAccessException e) {
                        throw new RuntimeException(e);
                    } catch (InvocationTargetException e) {
                        throw new RuntimeException(e);
                    }
                }
            }
            out.println();
        }
    }

    private static class MethodAndDescription implements Comparable<MethodAndDescription> {
        Method method;
        MDCDescription description;

        public MethodAndDescription(Method method, MDCDescription description) {
            this.method = method;
            this.description = description;
        }

        @Override
        public int compareTo(MethodAndDescription o) {
            return this.description.pos() - o.description.pos();
        }

        public MDCDescription getDescription() {
            return description;
        }

        public Method getMethod() {
            return method;
        }
    }
}
