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

/**
 * Identifies a configuration parameter of a {@link Project}. A key has a name, a
 * {@link ConfigParamInfo} describing its value type, range and default, and may optionally be
 * restricted to certain goals.
 */
public interface ConfigKey {
	/**
	 * Returns the name of this configuration parameter.
	 *
	 * @return the parameter name
	 */
	public String getName();

	/**
	 * Returns the descriptor of this parameter's value type, valid range and default value.
	 *
	 * @return the parameter's value type descriptor
	 */
	public ConfigParamInfo<?> getInfo();

	/**
	 * Returns whether this parameter is applicable to the given goal; {@code true} by default.
	 *
	 * @param forGoal the goal to test applicability for
	 * @return {@code true} if this parameter applies to the given goal
	 */
	default public boolean isForGoal(GoalKey forGoal) {
		return true;
	}
}
