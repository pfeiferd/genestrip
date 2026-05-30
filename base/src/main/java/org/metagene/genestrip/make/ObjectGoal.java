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
package org.metagene.genestrip.make;

public abstract class ObjectGoal<T, P extends Project> extends Goal<P> {
	private T object;

	@SafeVarargs
	public ObjectGoal(P project, GoalKey key, Goal<P>... dependencies) {
		super(project, key, dependencies);
	}
	
	@Override
	protected void allDependentsMade() {
		cleanThis();
	}
	
	public final synchronized T get() {
		make();
		return object;
	}

	protected final void set(T object) {
		if (getLogger().isTraceEnabled()) {
			getLogger().trace("Setting " + this + " to " + (object == null ? null : logSetObject(object)));
		}
		this.object = object;
	}
	
	protected String logSetObject(T object) {
		return object.toString();
	}

	@Override
	protected void doCleanThis() {
		set(null);
	}

	@Override
	public final boolean isMade() {
		return object != null;
	}

	@Override
	public void dump() {
		doCleanThis();
	}
	
	@Override
	public String toString() {
		return "object goal: " + getKey().getName();
	}
}
