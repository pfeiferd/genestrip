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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.logging.Log;
import org.metagene.genestrip.util.GSLogFactory;

public abstract class Goal<P extends Project> {
	private final Log logger;
	private final GoalKey goalKey;
	private final Goal<P>[] dependencies;
	private final P project;
	private final List<Goal<P>> dependents;

	@SafeVarargs
	public Goal(P project, GoalKey goalKey, Goal<P>... dependencies) {
		this.goalKey = initKey(goalKey);
		logger = GSLogFactory.getLog(this.goalKey.getName());
		this.dependencies = dependencies;
		this.project = project;
		this.dependents = new ArrayList<>();
		for (Goal<P> dep : dependencies) {
			dep.addDependent(this);
		}
	}

	protected void addDependent(Goal<P> goal) {
		dependents.add(goal);
	}

	protected void dependentMade(Goal<P> goal) {
		boolean allMade = true;
		for (Goal<P> dep : dependents) {
			if (!dep.isMade()) {
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
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Making " + this);
				getLogger().info("Dependencies " + Arrays.toString(dependencies));
			}
			startMake();
			GSLogFactory.incN();
			for (Goal<P> dep : dependencies) {
				if (dep != null && !isWeakDependency(dep)) {
					dep.make();
				}
			}
			GSLogFactory.decN();
			// It is important to check "isMade()" again here because sometime the goal get
			// automatically
			// made via a dependent goal as there are "shortcuts".
			makeThis();
			endMake();
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Made " + this);
			}
		} else {
			alreadyMade();
		}
	}

	private void logHeapInfo() {
		if (getLogger().isTraceEnabled()) {
			long total = Runtime.getRuntime().totalMemory();
			long free = Runtime.getRuntime().freeMemory();
			getLogger().trace("Total heap size: " + (total / 1024 / 1024) + " MB");
			getLogger().trace("Used heap size: " + ((total - free) / 1024 / 1024) + " MB");
		}
	}

	protected void endMake() {
	}

	protected void startMake() {
	}

	protected void alreadyMade() {
	}

	public final void makeThis() {
		if (!isMade()) {
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Making this " + this);
			}
			logHeapInfo();
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
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Cleaning this " + this);
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
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Cleaning " + this);
			getLogger().info("Dependencies " + Arrays.toString(dependencies));
		}
		for (Goal<P> dep : dependencies) {
			if (enforceTransitive) {
				dep.clean(true);
			} else {
				dep.transitiveClean();
			}
		}
		cleanThis();
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Cleaned " + this);
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
}
