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

public abstract class Maker<P extends Project> {
	private final Log logger = GSLogFactory.getLog("maker");
	
	// Linked Map so as to be able to returns the goals in entry order...
	private final Map<GoalKey, Goal<P>> goalsByKey;
	private final P project;
	private GoalKey defaultGoal;

	public Maker(P project) {
		goalsByKey = new LinkedHashMap<GoalKey, Goal<P>>();
		this.project = project;
		createGoals();
	}
	
	protected P getProject() {
		return project;
	}

	protected Log getLogger() {
		return logger;
	}
	
	protected void registerDefaultGoal(Goal<P> goal) {
		registerGoal(goal);
		if (defaultGoal != null) {
			throw new IllegalStateException("duplicate default goal");			
		}
		defaultGoal = goal.getKey();
	}
	
	protected void registerGoal(Goal<P> goal) {
		if (goalsByKey.get(goal.getKey()) != null) {
			throw new IllegalStateException("duplicate goal name");
		}
		goalsByKey.put(goal.getKey(), goal);
	}
	
	public Set<GoalKey> getGoalKeys() {
		return goalsByKey.keySet();
	}
	
	// Returns goals in entry order...
	public Collection<Goal<P>> getGoals() {
		return goalsByKey.values();
	}
		
	public Goal<P> getGoal(String keyName) {
		return getGoal(getKeyByName(keyName));
	}
	
	public GoalKey getKeyByName(String keyName) {
		for (GoalKey key : getGoalKeys()) {
			if (key.getName().equals(keyName)) {
				return key;
			}
		}
		return null;
	}
	
	public Goal<P> getGoal(GoalKey key) {
		return goalsByKey.get(key);
	}
	
	public void make(GoalKey goalKey) {
		goalsByKey.get(goalKey).make();
	}

	public void cleanAll(GoalKey goalkey) {
		goalsByKey.get(goalkey).clean();
	}
	
	public void cleanTotal(GoalKey goalkey) {
		goalsByKey.get(goalkey).clean(true);
	}
	
	public void clean(GoalKey goalkey) {
		goalsByKey.get(goalkey).cleanThis();
	}
	
	public GoalKey getDefaultGoalKey() {
		return defaultGoal;
	}
	
	public void dump() {
		for (Goal<P> g : goalsByKey.values()) {
			g.dump();
		}
	}
	
	protected abstract void createGoals();
}