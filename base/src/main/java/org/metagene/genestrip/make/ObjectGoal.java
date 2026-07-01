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

/**
 * A {@link Goal} that produces an in-memory value of type {@code T} rather than files. The value is
 * computed lazily on the first {@link #get()} (which makes the goal) and held until the goal is
 * cleaned or dumped. Being treated as a "weak" dependency, an object goal is not made automatically
 * when a dependent is made and is cleaned once all its dependents have been made, so its value does
 * not linger in memory longer than needed.
 *
 * @param <T> the type of value produced by this goal
 * @param <P> the type of {@link Project} this goal belongs to
 */
public abstract class ObjectGoal<T, P extends Project> extends Goal<P> {
	private T object;

	/**
	 * Creates the object goal.
	 *
	 * @param project      the project this goal belongs to
	 * @param key          the key identifying this goal
	 * @param dependencies the goals this goal depends on
	 */
	@SafeVarargs
	public ObjectGoal(P project, GoalKey key, Goal<P>... dependencies) {
		super(project, key, dependencies);
	}
	
	/**
	 * Frees the held value once all dependents have been made, by cleaning this goal.
	 */
	@Override
	protected void allDependentsMade() {
		cleanThis();
	}
	
	/**
	 * Makes this goal if necessary and returns its computed value.
	 *
	 * @return the value produced by this goal
	 */
	public final synchronized T get() {
		make();
		return object;
	}

	/**
	 * Stores the given value as this goal's computed value.
	 *
	 * @param object the value to store, may be {@code null}
	 */
	protected final void set(T object) {
		if (getLogger().isTraceEnabled()) {
			getLogger().trace("Setting " + this + " to " + (object == null ? null : logSetObject(object)));
		}
		this.object = object;
	}
	
	/**
	 * Returns the string used to represent the value in trace logging when it is set; may be
	 * overridden to avoid logging large objects.
	 *
	 * @param object the value being set
	 * @return the string representation used for trace logging
	 */
	protected String logSetObject(T object) {
		return object.toString();
	}

	@Override
	protected void doCleanThis() {
		set(null);
	}

	/**
	 * Returns whether the value has been computed, i.e. is non-{@code null}.
	 */
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
