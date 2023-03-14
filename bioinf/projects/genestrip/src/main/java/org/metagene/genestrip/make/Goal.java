package org.metagene.genestrip.make;

import java.util.Arrays;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

public abstract class Goal {
	private final Log logger = LogFactory.getLog(getClass());

	private final String name;
	private final Goal[] dependencies;

	public Goal(String name, Goal... dependencies) {
		this.name = name;
		this.dependencies = dependencies;
	}

	public String getName() {
		return name;
	}

	public Goal[] getDependencies() {
		return dependencies;
	}

	protected Log getLogger() {
		return logger;
	}

	public abstract boolean isMade();

	public final void make() {
		if (!isMade()) {
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Making " + this);
				getLogger().info("Dependencies " + Arrays.toString(dependencies));
			}
			startMake();
			for (Goal dep : dependencies) {
				if (!dep.isMade()) {
					dep.make();
				}
			}
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Making this " + this);
			}
			makeThis();
			endMake();
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Made " + this);
			}
		}
	}

	protected void endMake() {
	}

	protected void startMake() {
	}

	public void makeThis() {
	}

	public void cleanThis() {
	}

	@Override
	public String toString() {
		return "goal: " + getName();
	}

	public final void clean() {
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Cleaning " + this);
			getLogger().info("Dependencies " + Arrays.toString(dependencies));
		}
		for (Goal dep : dependencies) {
			dep.clean();
		}
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Cleaning this " + this);
		}
		cleanThis();
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Cleaned " + this);
		}
	}
}
