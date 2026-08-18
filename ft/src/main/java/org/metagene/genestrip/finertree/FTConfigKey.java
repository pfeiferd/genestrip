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
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.make.MDDescription;
import org.metagene.genestrip.make.ConfigKey;
import org.metagene.genestrip.make.ConfigParamInfo;
import org.metagene.genestrip.make.GoalKey;
import org.metagene.genestrip.probfilter.BlockedBloomFilter;
import org.metagene.genestrip.tax.Rank;
import org.metagene.genestrip.tax.TaxIdInfo;

import java.io.PrintStream;
import java.lang.annotation.Annotation;
import java.util.*;

/**
 * Enumeration of the configuration parameters added by the finer-tree (FT) extension. Each constant
 * binds a name to its {@link ConfigParamInfo} and the goals it applies to.
 */
public enum FTConfigKey implements ConfigKey {
    /** The cluster distance method used when performing agglomerative clustering. */
    @MDDescription("The cluster distance to be used when performing agglomerative clustering.")
    CLUSTER_METHOD("clusterMethod", new MethodConfigParamInfo(SimpleAggloClustering.Method.SINGLE_LINKAGE), FTGoalKey.DENDROGRAM),
    /** Whether the clustering merges on the Jaccard index rather than on containment. */
    @MDDescription("Whether the similarity of two of a node's children is their Jaccard index - the *k*-mers they share over the *k*-mers either has - or, with `false`, their containment: the shared *k*-mers over the smaller of the two sets. "
            + "It is what the clustering of `dendrogram` merges on. "
            + "**The rule: `false` where the children are individual genomes of comparable true size but differing assembly quality, `true` where they are taxa.** "
            + "A draft assembly's *k*-mer set is small because sequence is missing rather than different, and Jaccard reads the two alike, so a dendrogram over drafts clusters by assembly completeness; containment is invariant to that. "
            + "Between a genus and its species the children are taxa whose *k*-mer sets legitimately differ by orders of magnitude, and containment would score a sparse child as identical to whichever dense one contains it.")
    JACCARD_SIM("jaccardSim", new ConfigParamInfo.BooleanConfigParamInfo(true), FTGoalKey.INTERSECT_COUNT),
    /** Whether the Jaccard-index denominator includes the k-mer counts of all descendants. */
    @MDDescription("Whether to include the *k*-mer counts of *all* descendents for any two considered species in the denominator of the Jaccard-index. If not, only the *k*-mer counts right for the two considered species are used.")
    WITH_DESCENDANT_COUNTS("withDescendantCounts", new ConfigParamInfo.BooleanConfigParamInfo(false), FTGoalKey.INTERSECT_COUNT),
    /** Whether the LaTeX dendrogram is turned so that species names align horizontally. */
    @MDDescription("Whether the dendrogram in LaTeX has the species names aligned horizontally (with the entired diagram turned) or not.")
    TURN_LATEX("turnLatex", new ConfigParamInfo.BooleanConfigParamInfo(true), FTGoalKey.DENDRO_LATEX),
    /** Stretch factor for the dendrogram in TikZ's native x coordinate. */
    @MDDescription("The factory for stretching the dendrogram in TikZ's native *x* coordinate.")
    X_FACTOR_LATEX("xFactorLatex", new ConfigParamInfo.DoubleConfigParamInfo(0, Double.MAX_VALUE, 1), FTGoalKey.DENDRO_LATEX),
    /** Stretch factor for the dendrogram in TikZ's native y coordinate. */
    @MDDescription("The factory for stretching the dendrogram in TikZ's native *y* coordinate.")
    Y_FACTOR_LATEX("yFactorLatex", new ConfigParamInfo.DoubleConfigParamInfo(0, Double.MAX_VALUE, 8), FTGoalKey.DENDRO_LATEX),
    /** Value of {@code scale} in the {@code tikzpicture} environment of a LaTeX dendrogram. */
    @MDDescription("The factor for `scale` in the 'tikzpicture' environment of a LaTex dendrogram.")
    TIKZ_SCALE_FACTOR("tikzScaleFactor", new ConfigParamInfo.DoubleConfigParamInfo(0, Double.MAX_VALUE, 1), FTGoalKey.DENDRO_LATEX),
    /** Whether dendrogram similarities are scaled logarithmically. */
    @MDDescription("Whether to do logarithmic scaling of the similarity `sim` in dendrograms (via `1 - log (sim) / log (min_sim)`.")
    SIM_LOG_SCALING("simLogScaling", new ConfigParamInfo.BooleanConfigParamInfo(false), FTGoalKey.DENDRO_LATEX),
    /** How the Bloom filter of {@code kmerindexbloom} is sized. */
    @MDDescription("How the Bloom filter of `kmerindexbloom` is sized. `upperBound` uses the conservative bound, "
            + "which assumes every *k*-mer of a taxon to occur in every one of its genomes and needs nothing read "
            + "beforehand; genomes of one taxon share most of their *k*-mers, so that bound can exceed the truth by "
            + "orders of magnitude and the filter with it. `distinct` uses the estimate of `kmerindexsize`, which "
            + "reads the sequences once more to sketch the pairs. `auto` uses that estimate and falls back to the "
            + "bound, reading a second time, should the estimate turn out to have fallen short.")
    FT_BLOOM_FILTER_SIZING("ftBloomFilterSizing", new GSConfigKey.BloomFilterSizingConfigParamInfo(GSConfigKey.BloomFilterSizing.AUTO), FTGoalKey.KMER_INDEX_BLOOM),
    /** False positive probability of the {@code kmerindexbloom} filter. */
    @MDDescription("False positive probability (FPP) of the *k*-mer index filter built by `kmerindexbloom`. "
            + "At the default and above, the filter is a `BlockedBloomFilter`, which reaches about 1.3 per cent "
            + "once it holds what it was sized for and is the faster of the two; below it, a plain `BloomFilter` "
            + "is used instead, since the blocked one sets a fixed four bits per key and so cannot reach a lower "
            + "rate without spending far more memory than the optimal sizing needs. "
            + "WHY IT MATTERS FOR REFINEMENT: `ftupdatedb` places a *k*-mer at the smallest refined node covering "
            + "every child the filter reports it in, so a single false positive pushes it up the tree. The filter "
            + "is asked once per direct child, which makes the rate that counts `1 - (1 - fpp)^(C + 1)` for a node "
            + "with `C` children: harmless at twenty, but at the default FPP a node with 314 children leaves only "
            + "1.6 per cent of its *k*-mers placeable at all. Pick the FPP from the largest child count a refined "
            + "node has, roughly `-ln(P) / C` for a share `P` of *k*-mers to stay unaffected; at `C = 315` and "
            + "`P = 0.95` that is about `1e-4`. Memory grows only with `log2(1 / fpp)`, so this costs about twice "
            + "the bits per entry, and note that the filter has to be rebuilt - delete `<db>_storekmerindex.ser.gz`, "
            + "which a stale run would otherwise load at the FPP it was written with.")
    FT_INDEX_BLOOM_FILTER_FPP("ftIndexBloomFilterFpp", new ConfigParamInfo.DoubleConfigParamInfo(0, 1, BlockedBloomFilter.DEFAULT_FPP, true), FTGoalKey.KMER_INDEX_BLOOM),
    /** One in how many k-mers {@code kmerindexsize} looks at, or 1 to look at all of them. */
    @MDDescription("One in how many *k*-mers the goal `kmerindexsize` looks at when it estimates how many "
            + "(*k*-mer, leaf) pairs the filter of `kmerindexbloom` will hold; `1` looks at all of them. The "
            + "*k*-mers are selected by the *k*-mer itself, so one is either always looked at or never, however "
            + "often it is met - which is what makes the count of the sample scale to the whole. Everything not "
            + "selected is skipped before the store is consulted, which is where that pass spends its time. "
            + "Measured against the true number of pairs of real RefSeq sequence, one in 16 was off by less than "
            + "0.2 per cent and one in 64 by half a per cent, against a Bloom filter sizing that tolerates a few "
            + "per cent. Note that a sampled estimate is no longer exact for a small index, where the sketch "
            + "would otherwise have counted precisely; `upperBound` sizing ignores this setting altogether.")
    FT_KMER_INDEX_SIZE_SAMPLING("ftKMerIndexSizeSampling", new ConfigParamInfo.IntConfigParamInfo(1, Integer.MAX_VALUE, 16), FTGoalKey.KMER_INDEX_SIZE),
    /** Maximum number of dendrograms put into one LaTeX file by the {@code allinonelatex} goal. */
    @MDDescription("Maximum number of dendrograms put in one LaTeX file via the goal `allinonelatex`.")
    ALLINONE_CHUNK_SIZE("allInOneChunkSize", new ConfigParamInfo.IntConfigParamInfo(1, Integer.MAX_VALUE, 50)),
    /** The ranks or tax ids at which the taxonomy tree is refined, as {@link RefinementPosition}s. */
    @MDDescription("The ranks or tax ids for which the taxonomy tree is supposed to be refined.")
    REFINEMENT_POSITIONS("refinementPositions", new ConfigParamInfo.ListConfigParamInfo<>(Collections.unmodifiableList(Collections.singletonList(RefinementPosition.ALL_POSITIONS)
            // This was the old default here:
            /* Arrays.asList(new RefinementPosition(Rank.GENUS), new RefinementPosition(Rank.SPECIES_GROUP), new RefinementPosition(Rank.SUBGENUS))*/)) {
        @Override
        public String getTypeDescriptor() {
            return "comma-separated list of values of `<rank>`, `<taxid>` or else `*` which means all ranks and taxids are included. `>`, `<`, `>=`, `<=` and `=` may precede a rank which means nodes above, below the given rank etc. are included";
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

    /**
     * Returns whether this configuration key is internal.
     *
     * @return whether this configuration key is internal and therefore hidden from user documentation
     */
    public boolean isInternal() {
        return internal;
    }

    /**
     * Returns the textual name of this configuration key as used in property files and on the
     * command line.
     *
     * @return the configuration key's name
     */
    @Override
    public String getName() {
        return name;
    }

    /**
     * Returns the parameter descriptor holding the type, value range and default value of this
     * configuration key.
     *
     * @return the {@link ConfigParamInfo} associated with this key
     */
    public ConfigParamInfo<?> getInfo() {
        return param;
    }

    /**
     * @return whether this configuration key applies to the given goal; a {@code null} goal or a key
     *         declared for no specific goal matches any goal
     */
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

    /**
     * Parameter descriptor for a {@link SimpleAggloClustering.Method}, i.e. the cluster-distance
     * method used in agglomerative clustering. Values are nominal and parsed from the method's enum
     * name.
     */
    public static class MethodConfigParamInfo extends ConfigParamInfo<SimpleAggloClustering.Method> {
        /**
         * Creates a method parameter descriptor with the given default clustering method.
         *
         * @param defaultValue the default clustering method
         */
        public MethodConfigParamInfo(SimpleAggloClustering.Method defaultValue) {
            super(defaultValue);
        }

        /**
         * Tests whether the given object is an acceptable value, i.e. {@code null} or a
         * {@link SimpleAggloClustering.Method}.
         *
         * @param o the value to check
         * @return whether the value is {@code null} or a clustering method
         */
        @Override
        public boolean isValueInRange(Object o) {
            return o == null || o instanceof SimpleAggloClustering.Method;
        }

        /**
         * Parses a clustering method from its enum name.
         *
         * @param s the method name to parse
         * @return the corresponding {@link SimpleAggloClustering.Method}
         */
        @Override
        protected SimpleAggloClustering.Method fromString(String s) {
            return SimpleAggloClustering.Method.valueOf(s);
        }

        /**
         * Returns a Markdown fragment listing all available clustering method names as the value
         * range of this parameter.
         *
         * @return a comma-separated Markdown list of the available clustering method names
         */
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

        /**
         * Returns the type descriptor of this parameter.
         *
         * @return the string {@code "nominal"}
         */
        @Override
        public String getTypeDescriptor() {
            return "nominal";
        }
    }

    /**
     * Prints a Markdown table describing all non-internal FT configuration keys that apply to the
     * given goal (or all such keys if {@code filterGoalKey} is {@code null}).
     *
     * @param ps            the stream the Markdown table is written to
     * @param filterGoalKey the goal to filter keys by, or {@code null} to include all keys
     */
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

    /**
     * Specifies a position in the taxonomy tree at which the finer tree should be refined, given
     * either as a taxonomic {@link Rank} together with a comparison {@link Limit} or as a specific
     * tax id. Parsed from tokens such as {@code >=species} or {@code *}.
     */
    public static class RefinementPosition {
        /**
         * Sentinel refinement position that matches every taxonomy node.
         */
        public static final RefinementPosition ALL_POSITIONS = new RefinementPosition();

        /**
         * The comparison operator applied to a rank when matching nodes, with its textual token.
         */
        public enum Limit {
            /** Equality comparison, token {@code =}. */
            EQ("="),
            /** Strictly-larger comparison, token {@code >}. */
            LARGER(">"),
            /** Strictly-less comparison, token {@code <}. */
            LESS("<"),
            /** Larger-or-equal comparison, token {@code >=}. */
            LARGER_EQ(">="),
            /** Less-or-equal comparison, token {@code <=}. */
            LESS_EQ("<="),
            /** Match-all comparison, token {@code *}. */
            ALL("*");
            private final String comp;

            Limit(String comp) {
                this.comp = comp;
            }

            /**
             * Returns the textual token of this comparison operator.
             *
             * @return the textual token of this comparison operator (e.g. {@code >=})
             */
            public String getComp() {
                return comp;
            }

            /**
             * Returns the comparison limit whose token the given string starts with.
             *
             * @param token the string to inspect
             * @return the limit whose token the given string starts with, or {@code null} if none
             */
            public static Limit fromString(String token) {
                // Pick the longest matching token so that two-character operators (">=", "<=") win
                // over their single-character prefixes (">", "<") regardless of enum order.
                Limit best = null;
                for (Limit limit : Limit.values()) {
                    if (token.startsWith(limit.comp)
                            && (best == null || limit.comp.length() > best.comp.length())) {
                        best = limit;
                    }
                }
                return best;
            }
        }

        private final Limit limit;
        private final String taxid;
        private final Rank rank;

        private RefinementPosition() {
            taxid = null;
            rank = null;
            limit = Limit.ALL;
        }

        /**
         * Creates a refinement position matching exactly the given rank.
         *
         * @param upperRank the rank to match at
         */
        public RefinementPosition(Rank upperRank) {
            this(upperRank, null, Limit.EQ);
        }

        /**
         * Creates a refinement position for exactly one of a rank or a tax id.
         *
         * @param rank  the rank to match at, or {@code null} if a tax id is given
         * @param taxid the tax id to match, or {@code null} if a rank is given
         * @param limit the comparison operator applied to the rank
         * @throws IllegalArgumentException if not exactly one of {@code rank}/{@code taxid} is given,
         *         if {@code limit} is {@code null}, if a tax id is combined with a limit other than
         *         {@link Limit#EQ}, or if the rank is indeterminate
         */
        public RefinementPosition(Rank rank, String taxid, Limit limit) {
            if ((taxid != null &&  rank != null) || (taxid == null && rank == null) ) {
                throw new IllegalArgumentException("exactly one argument must not be null");
            }
            if (limit == null) {
                throw new IllegalArgumentException("limit must not be null");
            }
            if (taxid != null && !Limit.EQ.equals(limit)) {
                throw new IllegalArgumentException("limit for taxid must be EQ");
            }
            if (rank != null && rank.isIndeterminate()) {
                throw new IllegalArgumentException("rank must not be indeterminate");
            }
            this.rank = rank;
            this.taxid = taxid;
            this.limit = limit;
        }

        /**
         * Returns the comparison operator of this position.
         *
         * @return the comparison operator of this position
         */
        public Limit getLimit() {
            return limit;
        }

        /**
         * Returns the rank of this position.
         *
         * @return the rank of this position, or {@code null} if it is defined by a tax id
         */
        public Rank getRank() {
            return rank;
        }

        /**
         * Returns the tax id of this position.
         *
         * @return the tax id of this position, or {@code null} if it is defined by a rank
         */
        public String getTaxid() {
            return taxid;
        }

        @Override
        public String toString() {
            String s = limit.getComp();
            return s + (taxid == null ? (rank == null ? "" : rank.toString()) : taxid);
        }

        /**
         * Tests whether the given taxonomy node matches this refinement position.
         *
         * @param node the taxonomy node to test
         * @return whether the given node matches this position according to its rank/tax id and
         *         comparison limit
         */
        public boolean isMatchingNodeForPosition(TaxIdInfo node) {
            switch (limit) {
                case ALL:
                    return true;
                case LARGER_EQ:
                    if (rank.ordinal() == node.getRankOrdinal()) {
                        return true;
                    }
                case LARGER:
                    if (node.getRank().isIndeterminate()) {
                        while (node != null && rank.ordinal() != node.getRankOrdinal()) {
                            node = node.getParent();
                        }
                        return node == null;
                    } else {
                        return rank.ordinal() > node.getRankOrdinal();
                    }
                case LESS_EQ:
                    if (rank.ordinal() == node.getRankOrdinal()) {
                        return true;
                    }
                case LESS:
                    if (node.getRank().isIndeterminate()) {
                        while (node != null && rank.ordinal() != node.getRankOrdinal()) {
                            node = node.getParent();
                        }
                        return node != null;
                    } else {
                        return rank.ordinal() < node.getRankOrdinal();
                    }
                default:
                    if (rank != null) {
                        return rank.ordinal() == node.getRankOrdinal();
                    } else {
                        return taxid.equals(node.getTaxId());
                    }
            }
        }

        /**
         * Parses a refinement position from a token, e.g. {@code *} (all positions), {@code >=species}
         * (a rank with a comparison limit) or a bare tax id.
         *
         * @param token the token to parse
         * @return the parsed refinement position
         */
        public static RefinementPosition valueOf(String token) {
            if (Limit.ALL.getComp().equals(token)) {
                return ALL_POSITIONS;
            }
            token = token.trim();
            Limit limit = Limit.fromString(token);
            if (limit != null) {
                token = token.substring(limit.getComp().length());
            }
            else {
                limit = Limit.EQ;
            }
            Rank rank = Rank.byName(token);
            String taxid = rank == null ? token : null;
            return new RefinementPosition(rank, taxid, limit);
        }

        /**
         * Finds the first refinement position in the collection that matches the given node.
         *
         * @param node                the taxonomy node to match
         * @param refinementPositions the refinement positions to search
         * @return the first refinement position in the collection that matches the given node, or
         *         {@code null} if none matches
         */
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
