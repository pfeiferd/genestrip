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

import java.util.*;

import org.apache.commons.logging.Log;
import org.metagene.genestrip.util.GSLogFactory;

/**
 * Base class of the make-like build framework: a goal produces some output (files or an in-memory
 * value) and declares dependencies on other goals. Making a goal first makes its (non-weak)
 * prerequisites and is skipped when the goal is already made and up to date. Subclasses define how
 * the output is produced ({@link #doMakeThis()}), removed ({@link #doCleanThis()}) and detected
 * ({@link #isMade()}).
 *
 * @param <P> the type of {@link Project} this goal belongs to
 */
public abstract class Goal<P extends Project> {
	private final Log logger;
	private final GoalKey goalKey;
	private final Goal<P>[] dependencies;
	private final P project;
	private final List<Goal<P>> dependents;

	private boolean potentiallyRequiredForMake;
	private boolean cleaning;

	/**
	 * Creates a goal for the given project and key, wiring up the given dependencies and registering
	 * this goal as a dependent of each of them.
	 *
	 * @param project      the project this goal belongs to
	 * @param goalKey      the key identifying this goal
	 * @param dependencies the goals this goal depends on
	 * @throws IllegalStateException if the dependencies would introduce a cyclic dependency
	 */
	@SafeVarargs
	public Goal(P project, GoalKey goalKey, Goal<P>... dependencies) {
		this.goalKey = initKey(goalKey);
		logger = GSLogFactory.getLog(this.goalKey.getName());
		this.dependencies = dependencies;
		this.project = project;
		if (hasTransDependencyFor(this)) {
			throw new IllegalStateException("Cyclic dependency found for: " + this);
		}
		this.dependents = new ArrayList<>();
		for (Goal<P> dep : dependencies) {
			if (dep != null) {
				dep.addDependent(this);
			}
		}
	}

	private void addDependent(Goal<P> goal) {
		dependents.add(goal);
	}

	private void removeDependent(Goal<P> goal) {
		dependents.remove(goal);
	}

	// Reverses the dependent registration done in the constructor. Intended for short-lived goals
	// (e.g. Maker's internal goals): detaching them keeps the dependent lists of long-lived shared
	// goals from growing across repeated make()/clean() calls and lets the transient goal be GC'd.
	/**
	 * Detaches this goal from its dependencies by reversing the dependent registration done in the
	 * constructor.
	 */
	protected void detachFromDependencies() {
		for (Goal<P> dep : dependencies) {
			if (dep != null) {
				dep.removeDependent(this);
			}
		}
	}

	/**
	 * Returns whether the given goal is a transitive (direct or indirect) dependency of this goal.
	 *
	 * @param candidate the goal to test
	 * @return {@code true} if the candidate is a transitive dependency of this goal
	 */
	public boolean hasTransDependencyFor(Goal<P> candidate) {
		for (Goal<P> dep : dependencies) {
			if (dep != null) {
				if (dep.equals(candidate) || dep.hasTransDependencyFor(candidate)) {
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Callback invoked on a dependency when one of its dependents has been made; once all
	 * potentially required dependents are made, {@link #allDependentsMade()} is triggered.
	 *
	 * @param goal the dependent goal that has just been made
	 */
	protected void dependentMade(Goal<P> goal) {
		boolean allMade = true;
		for (Goal<P> dep : dependents) {
			if (dep.isPotentiallyRequired() && !dep.isMade()) {
				allMade = false;
				break;
			}
		}
		if (allMade) {
			allDependentsMade();
		}
	}

	/**
	 * Hook invoked once all potentially required dependents of this goal have been made; empty by
	 * default and meant to be overridden.
	 */
	// For override:
	protected void allDependentsMade() {
	}

	// Chance to change the key by overriding this method.
	/**
	 * Returns the key to actually use for this goal, giving subclasses a chance to change it by
	 * overriding this method.
	 *
	 * @param keyFromConstructor the key passed to the constructor
	 * @return the key to use for this goal
	 */
	protected GoalKey initKey(GoalKey keyFromConstructor) {
		return keyFromConstructor;
	}

	/**
	 * Returns the project this goal belongs to.
	 *
	 * @return the project this goal belongs to
	 */
	public P getProject() {
		return project;
	}

	/**
	 * Returns the key identifying this goal.
	 *
	 * @return the key identifying this goal
	 */
	public GoalKey getKey() {
		return goalKey;
	}

	/**
	 * Returns the goals this goal directly depends on.
	 *
	 * @return the goals this goal directly depends on
	 */
	public Goal<P>[] getDependencies() {
		return dependencies;
	}

	/**
	 * Returns the logger for this goal.
	 *
	 * @return the logger for this goal
	 */
	protected Log getLogger() {
		return logger;
	}

	/**
	 * Returns whether the given dependency is "weak", i.e. an {@link ObjectGoal}, and therefore not
	 * automatically made when this goal is made.
	 *
	 * @param toGoal the dependency goal to test
	 * @return {@code true} if the dependency is weak and not made automatically
	 */
	public boolean isWeakDependency(Goal<P> toGoal) {
		return toGoal instanceof ObjectGoal;
	}

	/**
	 * Returns whether this goal's output already exists and is up to date, so that making it can be
	 * skipped.
	 *
	 * @return {@code true} if this goal is already made and up to date
	 */
	public abstract boolean isMade();

	/**
	 * Makes this goal: if it is not already made, marks it (and its dependencies) as potentially
	 * required, makes all non-weak dependencies, then makes this goal itself. Does nothing if the
	 * goal is already made.
	 */
	public final void make() {
		if (!isMade()) {
			if (getLogger().isDebugEnabled()) {
				getLogger().debug("Making " + this);
				getLogger().debug("Dependencies " + Arrays.toString(dependencies));
			}
			try {
				markPotentiallyRequired(true);
				startMake();
				GSLogFactory.incN();
				for (Goal<P> dep : dependencies) {
					if (dep != null && !isWeakDependency(dep)) {
						dep.make();
					}
				}
				GSLogFactory.decN();
				makeThis();
				endMake();
			} finally {
				markPotentiallyRequired(false);
			}
			if (getLogger().isDebugEnabled()) {
				getLogger().debug("Made " + this);
			}
		} else {
			alreadyMade();
		}
	}

	/**
	 * Marks or unmarks this goal and all its dependencies as potentially required for making.
	 *
	 * @param mark {@code true} to mark as potentially required, {@code false} to unmark
	 */
	protected void markPotentiallyRequired(boolean mark) {
		// An simple optimization: We don't need to do things another time -
		// neither on this goal nor on its dependencies.
		if (potentiallyRequiredForMake == mark) {
			return;
		}
		potentiallyRequiredForMake = mark;
		for (Goal<P> dep : dependencies) {
			if (dep != null) {
				dep.markPotentiallyRequired(mark);
			}
		}
	}

	/**
	 * Returns whether this goal is currently marked as potentially required for making.
	 *
	 * @return {@code true} if this goal is marked as potentially required
	 */
	public boolean isPotentiallyRequired() {
		return potentiallyRequiredForMake;
	}

	/**
	 * Logs the current heap usage, at info level for {@link LogHeapInfo} goals and at debug level
	 * otherwise.
	 */
	protected void logHeapInfo() {
		if ((this instanceof LogHeapInfo && getLogger().isInfoEnabled())) {
			long total = Runtime.getRuntime().totalMemory();
			long free = Runtime.getRuntime().freeMemory();
			getLogger().info("Total heap size: " + (total / 1024 / 1024) + " MB");
			getLogger().info("Used heap size: " + ((total - free) / 1024 / 1024) + " MB");
		}
		else if (getLogger().isDebugEnabled()) {
			long total = Runtime.getRuntime().totalMemory();
			long free = Runtime.getRuntime().freeMemory();
			getLogger().debug("Total heap size: " + (total / 1024 / 1024) + " MB");
			getLogger().debug("Used heap size: " + ((total - free) / 1024 / 1024) + " MB");
		}
	}

	/**
	 * Hook invoked at the end of {@link #make()}, after this goal and its dependencies are made;
	 * empty by default.
	 */
	protected void endMake() {
	}

	/**
	 * Hook invoked at the start of {@link #make()}, before dependencies are made; empty by default.
	 */
	protected void startMake() {
	}

	/**
	 * Hook invoked by {@link #make()} when the goal is found to be already made; empty by default.
	 */
	protected void alreadyMade() {
	}

	/**
	 * Makes this goal's own output without making its dependencies, unless it is already made;
	 * notifies the dependencies afterwards via {@link #dependentMade(Goal)}.
	 */
	public final void makeThis() {
		// It is important to check "isMade()" again here because sometimes the goal gets
		// made automatically via a dependent goal as there are "shortcuts".
		if (!isMade()) {
			if (getLogger().isDebugEnabled()) {
				getLogger().debug("Making this " + this);
			}
			// logHeapInfo();
			doMakeThis();
			logHeapInfo();
			for (Goal<P> dep : dependencies) {
				if (dep != null) {
					dep.dependentMade(this);
				}
			}
		}
	}

	/**
	 * Performs the actual work of producing this goal's output; empty by default and overridden by
	 * subclasses.
	 */
	protected void doMakeThis() {
	}

	/**
	 * Cleans this goal's own output without cleaning its dependencies, unless it is already cleaned.
	 */
	public final void cleanThis() {
		if (!isCleaned()) {
			if (getLogger().isDebugEnabled()) {
				getLogger().debug("Cleaning this " + this);
			}
			doCleanThis();
		}
	}

	/**
	 * Performs the actual removal of this goal's output; empty by default and overridden by
	 * subclasses.
	 */
	protected void doCleanThis() {
	}

	/**
	 * Returns whether this goal is cleaned, i.e. its output no longer exists; by default the negation
	 * of {@link #isMade()}.
	 *
	 * @return {@code true} if this goal is cleaned (its output no longer exists)
	 */
	public boolean isCleaned() {
		return !isMade();
	}

	@Override
	public String toString() {
		return "goal: " + goalKey.getName();
	}

	/**
	 * Cleans this goal only if it permits transitive cleaning (see {@link #isAllowTransitiveClean()}).
	 */
	protected final void transitiveClean() {
		if (isAllowTransitiveClean()) {
			clean();
		}
	}

	/**
	 * Returns whether this goal may be cleaned transitively when a dependent goal is cleaned.
	 *
	 * @return {@code true} if this goal may be cleaned transitively
	 */
	public boolean isAllowTransitiveClean() {
		return goalKey.isTransClean();
	}

	/**
	 * Hook to release resources or detach transient state held by this goal; a no-op by default.
	 */
	public void dump() {
	}

	/**
	 * Cleans this goal and transitively cleans its dependencies that allow it.
	 */
	public final void clean() {
		clean(false);
	}

	/**
	 * Cleans this goal and then its dependencies. If {@code enforceTransitive} is {@code true} the
	 * dependencies are cleaned unconditionally, otherwise only those that allow transitive cleaning.
	 *
	 * @param enforceTransitive if {@code true}, dependencies are cleaned unconditionally; otherwise
	 *                          only those that allow transitive cleaning
	 */
	public final void clean(boolean enforceTransitive) {
		if (getLogger().isDebugEnabled()) {
			getLogger().debug("Cleaning " + this);
			getLogger().debug("Dependencies " + Arrays.toString(dependencies));
		}
		// Clean this first (cause the check if cleaned might depend on stuff
		// to be cleaned after).
		cleanThis();
		// Now clean dependencies.
		for (Goal<P> dep : dependencies) {
			if (dep != null) {
				if (enforceTransitive) {
					dep.clean(true);
				} else {
					dep.transitiveClean();
				}
			}
		}
		if (getLogger().isDebugEnabled()) {
			getLogger().debug("Cleaned " + this);
		}
	}

	/**
	 * Returns a new array containing the elements of {@code array} followed by the given
	 * {@code values}.
	 *
	 * @param <T>    the component type of the array
	 * @param array  the array whose elements come first
	 * @param values the values to append
	 * @return a new array containing the elements of {@code array} followed by {@code values}
	 */
	@SafeVarargs
	public static <T> T[] append(T[] array, T... values) {
		T[] result = Arrays.copyOf(array, array.length + values.length);
		for (int i = array.length; i < result.length; i++) {
			result[i] = values[i - array.length];
		}
		return result;
	}

	/**
	 * Returns the configuration value for the given key from this goal's project.
	 *
	 * @param key the configuration key
	 * @return the configuration value for the given key
	 */
	protected Object configValue(ConfigKey key) {
		return project.configValue(key);
	}

	/**
	 * Returns the configuration value for the given key as an {@code int}.
	 *
	 * @param key the configuration key
	 * @return the configuration value for the given key as an {@code int}
	 */
	protected int intConfigValue(ConfigKey key) {
		return project.intConfigValue(key);
	}

	/**
	 * Returns the configuration value for the given key as a {@code long}.
	 *
	 * @param key the configuration key
	 * @return the configuration value for the given key as a {@code long}
	 */
	protected long longConfigValue(ConfigKey key) {
		return project.longConfigValue(key);
	}

	/**
	 * Returns the configuration value for the given key as a {@code boolean}.
	 *
	 * @param key the configuration key
	 * @return the configuration value for the given key as a {@code boolean}
	 */
	protected boolean booleanConfigValue(ConfigKey key) {
		return project.booleanConfigValue(key);
	}

	/**
	 * Returns the configuration value for the given key as a {@code double}.
	 *
	 * @param key the configuration key
	 * @return the configuration value for the given key as a {@code double}
	 */
	protected double doubleConfigValue(ConfigKey key) {
		return project.doubleConfigValue(key);
	}

	/**
	 * Returns the configuration value for the given key as a {@code String}.
	 *
	 * @param key the configuration key
	 * @return the configuration value for the given key as a {@code String}
	 */
	protected String stringConfigValue(ConfigKey key) {
		return project.stringConfigValue(key);
	}

	// A marker interface regarding logging of memory on info level
	/**
	 * Marker interface for goals whose heap usage should be logged at info level.
	 */
	public interface LogHeapInfo {}
}
