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

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

public abstract class Maker<P> {
	private final Log logger = LogFactory.getLog(Maker.class);
	
	private final Map<String, Goal<P>> goalsByName;
	private String defaultGoalName;

	public Maker(P project) {
		goalsByName = new HashMap<String, Goal<P>>();
		createGoals(project);
	}

	protected Log getLogger() {
		return logger;
	}
	
	protected void registerDefaultGoal(Goal<P> goal) {
		registerGoal(goal);
		if (defaultGoalName != null) {
			throw new IllegalStateException("duplicate default goal");			
		}
		defaultGoalName = goal.getName();
	}
	
	protected void registerGoal(Goal<P> goal) {
		if (goalsByName.get(goal.getName()) != null) {
			throw new IllegalStateException("duplicate goal name");
		}
		goalsByName.put(goal.getName(), goal);
	}
	
	public Set<String> getGoalNames() {
		return goalsByName.keySet();
	}
	
	public Goal<P> getGoal(String name) {
		return goalsByName.get(name);
	}
	
	public void make(String goalName) {
		goalsByName.get(goalName).make();
	}

	public void cleanAll(String goalName) {
		goalsByName.get(goalName).clean();
	}
	
	public void cleanTotal(String goalName) {
		goalsByName.get(goalName).clean(true);
	}
	
	public void clean(String goalName) {
		goalsByName.get(goalName).cleanThis();
	}
	
	public String getDefaultGoalName() {
		return defaultGoalName;
	}
	
	protected abstract void createGoals(P project);
}