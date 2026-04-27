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
package org.metagene.genestrip.finertree.cluster;

public class DendrogramNode {
    public interface Visitor {
        public void preNode(DendrogramNode node);
        public void postNode(DendrogramNode node);
    }

    private final int valueIndex;
    private final DendrogramNode child1;
    private final DendrogramNode child2;
    private double similarity;
    private Object value;

    public DendrogramNode(DendrogramNode child1, DendrogramNode child2, double similarity) {
        this.valueIndex = -1;
        this.child1 = child1;
        this.child2 = child2;
        this.similarity = similarity;
    }

    public DendrogramNode(int valueIndex, double similarity) {
        this.valueIndex = valueIndex;
        this.child1 = null;
        this.child2 = null;
        this.similarity = similarity;
    }

    public int getValueIndex() {
        return valueIndex;
    }

    public DendrogramNode getChild1() {
        return child1;
    }

    public DendrogramNode getChild2() {
        return child2;
    }

    public double getSimilarity() {
        return similarity;
    }

    public void setValue(Object value) {
        this.value = value;
    }

    public Object getValue() {
        return value;
    }

    public void visit(Visitor visitor) {
        visitor.preNode(this);
        if (valueIndex == - 1) {
            child1.visit(visitor);
            child2.visit(visitor);
        }
        visitor.postNode(this);
    }

    public int size() {
        if (valueIndex == - 1) {
            return 1 + child1.size() + child2.size();
        }
        else {
            return 1;
        }
    }
}
