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
    @MDDescription("False positive probability of the Bloom filter used in Genestrip-FT.")
    FT_BLOOM_FILTER_FPP("ftBloomFilterFpp", new ConfigParamInfo.DoubleConfigParamInfo(0, 1, 0.001d), false, FTGoalKey.DB_QUALITY_COUNTS);

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

    public static class RefinementPosition {
        public static final RefinementPosition ALL_POSITIONS = new RefinementPosition();

        public enum Limit { EQ("="), LARGER(">"), LESS("<"), LARGER_EQ(">="), LESS_EQ("<="), ALL("*");
            private final String comp;

            Limit(String comp) {
                this.comp = comp;
            }

            public String getComp() {
                return comp;
            }

            public static Limit fromString(String token) {
                for (Limit limit : Limit.values()) {
                    if (token.startsWith(limit.comp)) {
                        return limit;
                    }
                }
                return null;
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

        public RefinementPosition(Rank upperRank) {
            this(upperRank, null, Limit.EQ);
        }

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

        public Limit getLimit() {
            return limit;
        }

        public Rank getRank() {
            return rank;
        }

        public String getTaxid() {
            return taxid;
        }

        @Override
        public String toString() {
            String s = limit.getComp();
            return s + (taxid == null ? (rank == null ? "" : rank.toString()) : taxid);
        }

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
