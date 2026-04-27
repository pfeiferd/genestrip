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
package org.metagene.genestrip.finertree;

import org.metagene.genestrip.finertree.cluster.SimpleAggloClustering;
import org.metagene.genestrip.make.MDDescription;
import org.metagene.genestrip.make.ConfigKey;
import org.metagene.genestrip.make.ConfigParamInfo;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxIdInfo;

import java.io.PrintStream;
import java.lang.annotation.Annotation;
import java.util.*;

public enum FTConfigKey implements ConfigKey {
    @MDDescription("The cluster distance to be used when performing agglomerative clustering.")
    CLUSTER_METHOD("clusterMethod", new MethodConfigParamInfo(SimpleAggloClustering.Method.SINGLE_LINKAGE), FTGoalKey.DENDROGRAM),
    @MDDescription("Whether to include the *k*-mer counts of *all* descendents for any two considered species in the denominator of the Jaccard-index. If not, only the *k*-mer counts right for the two considered species are used.")
    WITH_DESCENDANT_COUNTS("withDescendantCounts", new ConfigParamInfo.BooleanConfigParamInfo(false), FTGoalKey.INTERSECT_COUNT),
    @MDDescription("Whether the dendrogram in LaTeX has the species names aligned horizontally (with the entired diagram turned) or not.")
    TURN_LATEX("turnLatex", new ConfigParamInfo.BooleanConfigParamInfo(true), FTGoalKey.DENDRO_LATEX),
    @MDDescription("The factory for stretching the dendrogram in TikZ's native *x* coordinate.")
    X_FACTOR_LATEX("xFactorLatex", new ConfigParamInfo.DoubleConfigParamInfo(0, Double.MAX_VALUE, 1), FTGoalKey.DENDRO_LATEX),
    @MDDescription("The factory for stretching the dendrogram in TikZ's native *y* coordinate.")
    Y_FACTOR_LATEX("yFactorLatex", new ConfigParamInfo.DoubleConfigParamInfo(0, Double.MAX_VALUE, 8), FTGoalKey.DENDRO_LATEX),
    @MDDescription("The factor for `scale` in the 'tikzpicture' environment of a LaTex dendrogram.")
    TIKZ_SCALE_FACTOR("tikzScaleFactor", new ConfigParamInfo.DoubleConfigParamInfo(0, Double.MAX_VALUE, 1), FTGoalKey.DENDRO_LATEX),
    @MDDescription("Whether to do logarithmic scaling of the similarity `sim` in dendrograms (via `1 - log (sim) / log (min_sim)`.")
    SIM_LOG_SCALING("simLogScaling", new ConfigParamInfo.BooleanConfigParamInfo(false), FTGoalKey.DENDRO_LATEX),
    @MDDescription("False positive probability of the Bloom filter used in Genestrip-FT.")
    FT_BLOOM_FILTER_FPP("ftBloomFilterFpp", new ConfigParamInfo.DoubleConfigParamInfo(0, 1, 0.001d), false, FTGoalKey.KMER_INDEX_BLOOM),
    @MDDescription("Maximum number of dendrograms put in one LaTeX file via the goal `allinonelatex`.")
    ALLINONE_CHUNK_SIZE("allInOneChunkSize", new ConfigParamInfo.IntConfigParamInfo(1, Integer.MAX_VALUE, 50)),
    @MDDescription("The ranks or tax ids for which the taxonomy tree is supposed to be refined.")
    REFINEMENT_POSITIONS("refinementPositions", new ConfigParamInfo.ListConfigParamInfo<RefinementPosition>(Collections.unmodifiableList(Arrays.asList(new RefinementPosition(Rank.GENUS), new RefinementPosition(Rank.SPECIES_GROUP), new RefinementPosition(Rank.SUBGENUS)))) {
        @Override
        public String getTypeDescriptor() {
            return "comma-separated list of values of `<rank>`, `<taxid>` or else `*` which means all ranks and taxids are included";
        }

        @Override
        protected List<RefinementPosition> fromString(String qs) {
            List<RefinementPosition> res = new ArrayList<>();
            if (qs != null) {
                StringTokenizer tokenizer = new StringTokenizer(qs, ",;");
                while (tokenizer.hasMoreTokens()) {
                    RefinementPosition r = RefinementPosition.valueOf(tokenizer.nextToken().trim());
                    if (r != null) {
                        res.add(r);
                    }
                }
            }
            return res;
        }

        @Override
        public String getMDRangeDescriptor() {
            StringBuilder sb = new StringBuilder();
            sb.append("<rank> as subset of ");
            boolean first = true;
            for (Rank e : Rank.values()) {
                if (!first) {
                    sb.append(", ");
                }
                first = false;
                sb.append('`');
                sb.append(e.getName());
                sb.append('`');
            }
            return sb.toString();
        }
    }, FTGoalKey.KMER_INDEX_BLOOM, FTGoalKey.DENDRO_LATEX, FTGoalKey.INTERSECT_COUNT, FTGoalKey.INTERSECT_CSV);

    private final String name;
    private final ConfigParamInfo<?> param;
    private final boolean internal;
    private final FTGoalKey[] forGoals;

    FTConfigKey(String name, ConfigParamInfo<?> param, FTGoalKey... forGoals) {
        this(name, param, false, forGoals);
    }

    FTConfigKey(String name, ConfigParamInfo<?> param, boolean internal, FTGoalKey... forGoals) {
        this.name = name;
        this.param = param;
        this.internal = internal;
        this.forGoals = forGoals;
    }

    public boolean isInternal() {
        return internal;
    }

    @Override
    public String getName() {
        return name;
    }

    public ConfigParamInfo<?> getInfo() {
        return param;
    }

    public boolean isForGoal(GoalKey forGoal) {
        if (forGoal == null) {
            return true;
        }
        for (GoalKey id : forGoals) {
            if (forGoal.equals(id)) {
                return true;
            }
        }
        return false;
    }

    @Override
    public String toString() {
        return getName();
    }

    public static class MethodConfigParamInfo extends ConfigParamInfo<SimpleAggloClustering.Method> {
        public MethodConfigParamInfo(SimpleAggloClustering.Method defaultValue) {
            super(defaultValue);
        }

        @Override
        public boolean isValueInRange(Object o) {
            return o == null || o instanceof SimpleAggloClustering.Method;
        }

        @Override
        protected SimpleAggloClustering.Method fromString(String s) {
            return SimpleAggloClustering.Method.valueOf(s);
        }

        @Override
        public String getMDRangeDescriptor() {
            StringBuilder builder = new StringBuilder();
            SimpleAggloClustering.Method[] methods = SimpleAggloClustering.Method.values();
            for (int i = 0; i < methods.length; i++) {
                if (i > 0) {
                    builder.append(", ");
                }
                builder.append('`');
                builder.append(methods[i].name());
                builder.append('`');
            }
            return builder.toString();
        }

        @Override
        public String getTypeDescriptor() {
            return "nominal";
        }
    }

    public static void printMDConfigParamInfo(PrintStream ps, GoalKey filterGoalKey) {
        ps.print('|');
        ps.print("Name");
        ps.print('|');
        ps.print("Type");
        ps.print('|');
        ps.print("Value Range");
        ps.print('|');
        ps.print("Default");
        ps.print('|');
        ps.print("Description");
        ps.print('|');
        ps.print("For Goals");
        ps.print('|');
        ps.println();

        ps.print('|');
        ps.print('-');
        ps.print('|');
        ps.print('-');
        ps.print('|');
        ps.print('-');
        ps.print('|');
        ps.print('-');
        ps.print('|');
        ps.print('-');
        ps.print('|');
        ps.print('-');
        ps.print('|');
        ps.println();

        for (FTConfigKey configKey : FTConfigKey.values()) {
            if (!configKey.isInternal() && configKey.isForGoal(filterGoalKey)) {
                ps.print('|');
                ps.print('`');
                ps.print(configKey.getName());
                ps.print('`');
                ps.print('|');
                ps.print(configKey.getInfo().getTypeDescriptor());
                ps.print('|');
                ps.print(configKey.getInfo().getMDRangeDescriptor());
                ps.print('|');
                ps.print('`');
                ps.print(configKey.getInfo().getMDDefaultValue());
                ps.print('`');
                ps.print('|');
                Annotation[] annotations;
                try {
                    annotations = FTConfigKey.class.getField(configKey.name()).getAnnotations();
                } catch (NoSuchFieldException e) {
                    throw new RuntimeException(e);
                } catch (SecurityException e) {
                    throw new RuntimeException(e);
                }
                for (Annotation annotation : annotations) {
                    if (annotation instanceof MDDescription) {
                        ps.print(((MDDescription) annotation).value());
                        break;
                    }
                }
                ps.print('|');
                if (configKey.forGoals.length == 0) {
                    ps.print("all");
                } else {
                    boolean first = true;
                    for (GoalKey key : configKey.forGoals) {
                        if (!first) {
                            ps.print(", ");
                        }
                        first = false;
                        ps.print('`');
                        ps.print(key.getName());
                        ps.print('`');
                    }
                }
                ps.print('|');
                ps.println();
            }
        }
    }

    public static class RefinementPosition {
        public static RefinementPosition WILDCARD = new RefinementPosition();

        private final String upperTaxid;
        private final Rank upperRank;

        private RefinementPosition() {
            upperTaxid = null;
            upperRank = null;
        }

        public RefinementPosition(Rank upperRank) {
            this(upperRank, null);
        }

        public RefinementPosition(String taxid) {
            this(null, taxid);
        }

        public RefinementPosition(Rank upperRank, String taxid) {
            if ((taxid != null &&  upperRank != null) || (taxid == null && upperRank == null) ) {
                throw new IllegalArgumentException("exactly one argument must not be null");
            }
            this.upperRank = upperRank;
            this.upperTaxid = taxid;
        }

        public Rank getUpperRank() {
            return upperRank;
        }

        public String getUpperTaxid() {
            return upperTaxid;
        }

        @Override
        public String toString() {
            return upperTaxid == null ? (upperRank == null ? "*" : upperRank.toString()) : upperTaxid;
        }

        public boolean isMatchingNodeForPosition(TaxIdInfo node) {
            if (this == WILDCARD) {
                return true;
            }
            if (upperRank != null) {
                return upperRank.ordinal() == node.getRankOrdinal();
            } else {
                return upperTaxid.equals(node.getTaxId());
            }
        }

        public static RefinementPosition valueOf(String token) {
            String upperRankStr = token.trim();
            if ("*".equals(upperRankStr)) {
                return WILDCARD;
            }
            Rank upperRank = Rank.valueOf(upperRankStr);
            String taxid = upperRankStr == null ? upperRankStr : null;
            return new RefinementPosition(upperRank, taxid);
        }

        public static RefinementPosition getMatchingNodeFor(TaxIdInfo node, Collection<RefinementPosition> refinementPositions) {
            for (RefinementPosition interval : refinementPositions) {
                if (interval.isMatchingNodeForPosition(node)) {
                    return interval;
                }
            }
            return null;
        }
    }
}
