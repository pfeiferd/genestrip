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

import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.metagene.genestrip.util.GSLogFactory;

/**
 * Holds the goal graph of a {@link Project}: goals are created in {@link #createGoals()} and
 * registered by their {@link GoalKey}. Making or cleaning a set of goals is delegated to the goal
 * graph, optionally through a transient internal goal that aggregates the requested goals so shared
 * dependencies are made only once.
 *
 * @param <P> the type of {@link Project} this maker builds goals for
 */
public abstract class Maker<P extends Project> {
	private static GoalKey INTERNAL_KEY = new GoalKey.DefaultGoalKey("internalgoal");

	private final Log logger = GSLogFactory.getLog("maker");
	
	// Linked Map so as to be able to returns the goals in entry order...
	private final Map<GoalKey, Goal<P>> goalsByKey;
	private final P project;
	private GoalKey defaultGoal;

	/**
	 * Creates a maker for the given project and immediately builds its goal graph via
	 * {@link #createGoals()}.
	 *
	 * @param project the project whose goals are built
	 */
	public Maker(P project) {
		goalsByKey = new LinkedHashMap<GoalKey, Goal<P>>();
		this.project = project;
		createGoals();
	}
	
	/**
	 * Returns the project this maker builds goals for.
	 *
	 * @return the project
	 */
	protected P getProject() {
		return project;
	}

	/**
	 * Returns the logger of this maker.
	 *
	 * @return the maker logger
	 */
	protected Log getLogger() {
		return logger;
	}
	
	/**
	 * Registers the given goal and marks it as this maker's default goal.
	 *
	 * @param goal the goal to register as the default goal
	 * @throws IllegalStateException if a default goal has already been registered
	 */
	protected void registerDefaultGoal(Goal<P> goal) {
		registerGoal(goal);
		if (defaultGoal != null) {
			throw new IllegalStateException("duplicate default goal");			
		}
		defaultGoal = goal.getKey();
	}
	
	/**
	 * Registers the given goal under its key.
	 *
	 * @param goal the goal to register
	 * @throws IllegalStateException if a goal with the same key is already registered
	 */
	protected void registerGoal(Goal<P> goal) {
		if (goalsByKey.get(goal.getKey()) != null) {
			throw new IllegalStateException("duplicate goal name");
		}
		goalsByKey.put(goal.getKey(), goal);
	}
	
	/**
	 * Returns the keys of all registered goals.
	 *
	 * @return the registered goal keys
	 */
	public Set<GoalKey> getGoalKeys() {
		return goalsByKey.keySet();
	}

	/**
	 * Returns all registered goals in registration order.
	 *
	 * @return the registered goals in registration order
	 */
	// Returns goals in entry order...
	public Collection<Goal<P>> getGoals() {
		return goalsByKey.values();
	}

	/**
	 * Returns the goals registered under the given keys, in the given order.
	 *
	 * @param keys the keys of the goals to return
	 * @return the goals registered under the given keys, in the given order
	 */
	public Goal<P>[] getGoals(GoalKey... keys) {
		Goal<P>[] goals = new Goal[keys.length];
		for (int i = 0; i < keys.length; i++) {
			goals[i] = getGoal(keys[i]);
		}
		return goals;
	}
		
	/**
	 * Returns the goal whose key has the given name, or {@code null} if there is none.
	 *
	 * @param keyName the name of the goal key
	 * @return the matching goal, or {@code null} if there is none
	 */
	public Goal<P> getGoal(String keyName) {
		return getGoal(getKeyByName(keyName));
	}
	
	/**
	 * Returns the registered goal key with the given name, or {@code null} if there is none.
	 *
	 * @param keyName the name of the goal key
	 * @return the matching goal key, or {@code null} if there is none
	 */
	public GoalKey getKeyByName(String keyName) {
		for (GoalKey key : getGoalKeys()) {
			if (key.getName().equals(keyName)) {
				return key;
			}
		}
		return null;
	}
	
	/**
	 * Returns the goal registered under the given key, or {@code null} if there is none.
	 *
	 * @param key the goal key
	 * @return the registered goal, or {@code null} if there is none
	 */
	public Goal<P> getGoal(GoalKey key) {
		return goalsByKey.get(key);
	}

	/**
	 * Makes the given goals together (non-isolated), so shared dependencies are made only once.
	 *
	 * @param keys the keys of the goals to make
	 */
	public void make(GoalKey... keys) {
		make(false, keys);
	}

	/**
	 * Makes the given goals. If {@code isolate} is {@code true} each goal is made independently;
	 * otherwise they are made together through a single internal goal so shared dependencies are
	 * made only once.
	 *
	 * @param isolate whether each goal is made independently
	 * @param keys the keys of the goals to make
	 */
	public void make(boolean isolate, GoalKey... keys) {
		if (isolate) {
			for (GoalKey key : keys) {
				Goal goal = getGoal(key);
				if (goal != null) {
					goal.make();
				}
			}
		}
		else {
			runInternal(g -> g.make(), keys);
		}
	}

	/**
	 * Creates a transient internal goal that aggregates the goals for the given keys as its
	 * dependencies.
	 *
	 * @param keys the keys of the goals to aggregate
	 * @return the transient aggregating goal
	 */
	// We need this goal for the functionality regarding "potentiallyForMakeRequired" in class Goal
	protected Goal<P> createInternalGoal(GoalKey... keys) {
		return new Goal<P>(getProject(), INTERNAL_KEY, getGoals(keys)) {
			@Override
			public boolean isMade() {
				return false;
			}

			// This internal goal is created fresh for every make()/clean() call; dropping it from
			// its dependencies' dependent lists on dump() prevents those long-lived goals from
			// accumulating dead dependents. (Registered goals' dump() stays a no-op, so Maker.dump()
			// does not disturb the real inter-goal dependency edges.)
			@Override
			public void dump() {
				detachFromDependencies();
			}
		};
	}

	// Runs the action on a fresh internal goal and always detaches it afterwards (via dump()).
	private void runInternal(java.util.function.Consumer<Goal<P>> action, GoalKey... keys) {
		Goal<P> internalGoal = createInternalGoal(keys);
		try {
			action.accept(internalGoal);
		} finally {
			internalGoal.dump();
		}
	}

	/**
	 * Cleans the given goals and, transitively, their dependencies that allow transitive cleaning.
	 *
	 * @param keys the keys of the goals to clean
	 */
	public void cleanAll(GoalKey... keys) {
		runInternal(g -> g.clean(), keys);
	}

	/**
	 * Cleans the given goals and all their dependencies, enforcing transitive cleaning.
	 *
	 * @param keys the keys of the goals to clean
	 */
	public void cleanTotal(GoalKey... keys) {
		runInternal(g -> g.clean(true), keys);
	}

	/**
	 * Cleans only the given goals themselves, without cleaning their dependencies.
	 *
	 * @param keys the keys of the goals to clean
	 */
	public void clean(GoalKey... keys) {
		runInternal(g -> g.cleanThis(), keys);
	}
	
	/**
	 * Returns the key of this maker's default goal, or {@code null} if none was registered.
	 *
	 * @return the default goal key, or {@code null} if none
	 */
	public GoalKey getDefaultGoalKey() {
		return defaultGoal;
	}
	
	/**
	 * Dumps all registered goals, releasing their transient state (see {@link Goal#dump()}).
	 */
	public void dump() {
		for (Goal<P> g : goalsByKey.values()) {
			g.dump();
		}
	}
	
	/**
	 * Creates and registers all goals of this maker; called from the constructor.
	 */
	protected abstract void createGoals();
}