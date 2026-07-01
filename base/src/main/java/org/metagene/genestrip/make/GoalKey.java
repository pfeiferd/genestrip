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

import java.util.Objects;

/**
 * Identifies a goal within a {@link Maker}. A key has a name used for lookup and declares whether
 * goals identified by it may be cleaned transitively.
 */
public interface GoalKey {
	/**
	 * Returns the name used to look up the goal identified by this key.
	 *
	 * @return the goal name
	 */
	public String getName();
	/**
	 * Returns whether goals identified by this key may be cleaned transitively when one of their
	 * dependents is cleaned.
	 *
	 * @return {@code true} if transitive cleaning is allowed
	 */
	public boolean isTransClean();
	
	/**
	 * A simple {@link GoalKey} identified solely by its name, with transitive cleaning disabled.
	 */
	public static class DefaultGoalKey implements GoalKey {
		private final String name;
		
		/**
		 * Creates a key identified by the given name.
		 *
		 * @param name the goal name
		 */
		public DefaultGoalKey(String name) {
			this.name = name;
		}
		
		@Override
		public String getName() {
			return name;
		}

		@Override
		public boolean isTransClean() {
			return false;
		}

		@Override
		public String toString() {
			return name;
		}

		// Identity by name so a DefaultGoalKey resolves the same goal in Maker.goalsByKey
		// regardless of which instance is used. Only equal to other DefaultGoalKeys, never to the
		// GSGoalKey enum constants that share the map.
		@Override
		public boolean equals(Object obj) {
			if (this == obj) {
				return true;
			}
			if (!(obj instanceof DefaultGoalKey)) {
				return false;
			}
			return Objects.equals(name, ((DefaultGoalKey) obj).name);
		}

		@Override
		public int hashCode() {
			return Objects.hashCode(name);
		}
	}
}
