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
	
	public void clean(String goalName) {
		goalsByName.get(goalName).cleanThis();
	}
	
	public String getDefaultGoalName() {
		return defaultGoalName;
	}
	
	protected abstract void createGoals(P project);
}