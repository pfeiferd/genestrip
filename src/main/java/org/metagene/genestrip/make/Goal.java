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

import java.util.Arrays;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

public abstract class Goal<P> {
	private final Log logger;
	private final String name;
	private final Goal<P>[] dependencies;
	private final P project;

	@SafeVarargs
	public Goal(P project, String name, Goal<P>... dependencies) {
		logger = LogFactory.getLog(name);
		this.name = name;
		this.dependencies = dependencies;
		this.project = project;
	}

	public P getProject() {
		return project;
	}

	public String getName() {
		return name;
	}

	public Goal<P>[] getDependencies() {
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
			for (Goal<P> dep : dependencies) {
				if (dep != null && !dep.isMade()) {
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

	protected boolean isThisCleaned() {
		return false;
	}

	@Override
	public String toString() {
		return "goal: " + getName();
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
			}
			else {
				dep.transitiveClean();
			}
		}
		if (!isThisCleaned()) {
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Cleaning this " + this);
			}
			cleanThis();
			if (getLogger().isInfoEnabled()) {
				getLogger().info("Cleaned " + this);
			}
		}
	}
}
