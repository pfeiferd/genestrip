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

public abstract class Goal<P extends Project> {
	private final Log logger;
	private final GoalKey goalKey;
	private final Goal<P>[] dependencies;
	private final P project;
	private final List<Goal<P>> dependents;

	private boolean potentiallyForMakeRequired;

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
			dep.addDependent(this);
		}
	}

	private void addDependent(Goal<P> goal) {
		dependents.add(goal);
	}

	public boolean hasTransDependencyFor(Goal<P> candidate) {
		for (Goal<P> dep : dependencies) {
			if (dep.equals(candidate) || dep.hasTransDependencyFor(this)) {
				return true;
			}
		}
		return false;
	}

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

	// For override:
	protected void allDependentsMade() {
	}

	// Chance to change the key by overriding this method.
	protected GoalKey initKey(GoalKey keyFromConstructor) {
		return keyFromConstructor;
	}

	public P getProject() {
		return project;
	}

	public GoalKey getKey() {
		return goalKey;
	}

	public Goal<P>[] getDependencies() {
		return dependencies;
	}

	protected Log getLogger() {
		return logger;
	}

	public boolean isWeakDependency(Goal<P> toGoal) {
		return toGoal instanceof ObjectGoal;
	}

	public abstract boolean isMade();

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

	protected void markPotentiallyRequired(boolean mark) {
		potentiallyForMakeRequired = mark;
		for (Goal<P> dep : dependencies) {
			dep.markPotentiallyRequired(mark);
		}
	}

	public boolean isPotentiallyRequired() {
		return potentiallyForMakeRequired;
	}

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

	protected void endMake() {
	}

	protected void startMake() {
	}

	protected void alreadyMade() {
	}

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
				dep.dependentMade(this);
			}
		}
	}

	protected void doMakeThis() {
	}

	public final void cleanThis() {
		if (!isCleaned()) {
			if (getLogger().isDebugEnabled()) {
				getLogger().debug("Cleaning this " + this);
			}
			doCleanThis();
		}
	}

	protected void doCleanThis() {
	}

	public boolean isCleaned() {
		return !isMade();
	}

	@Override
	public String toString() {
		return "goal: " + goalKey.getName();
	}

	protected final void transitiveClean() {
		if (isAllowTransitiveClean()) {
			clean();
		}
	}

	public boolean isAllowTransitiveClean() {
		return true;
	}

	public void dump() {
	}
	
	public final void clean() {
		clean(false);
	}

	public final void clean(boolean enforceTransitive) {
		if (getLogger().isDebugEnabled()) {
			getLogger().debug("Cleaning " + this);
			getLogger().debug("Dependencies " + Arrays.toString(dependencies));
		}
		for (Goal<P> dep : dependencies) {
			if (enforceTransitive) {
				dep.clean(true);
			} else {
				dep.transitiveClean();
			}
		}
		cleanThis();
		if (getLogger().isDebugEnabled()) {
			getLogger().debug("Cleaned " + this);
		}
	}

	@SafeVarargs
	public static <T> T[] append(T[] array, T... values) {
		T[] result = Arrays.copyOf(array, array.length + values.length);
		for (int i = array.length; i < result.length; i++) {
			result[i] = values[i - array.length];
		}
		return result;
	}

	protected Object configValue(ConfigKey key) {
		return project.configValue(key);
	}

	protected int intConfigValue(ConfigKey key) {
		return project.intConfigValue(key);
	}

	protected long longConfigValue(ConfigKey key) {
		return project.longConfigValue(key);
	}

	protected boolean booleanConfigValue(ConfigKey key) {
		return project.booleanConfigValue(key);
	}

	protected double doubleConfigValue(ConfigKey key) {
		return project.doubleConfigValue(key);
	}

	protected String stringConfigValue(ConfigKey key) {
		return project.stringConfigValue(key);
	}

	// A marker interface regarding logging of memory on info level
	public interface LogHeapInfo {}
}
