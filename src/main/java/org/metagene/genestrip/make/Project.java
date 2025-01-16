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
import java.util.Properties;

import org.apache.commons.logging.Log;
import org.metagene.genestrip.GSConfigKey;
import org.metagene.genestrip.util.GSLogFactory;

public abstract class Project {
	private final Log logger;
	private final String name;
	private final Map<ConfigKey, Object> paramMap;
	
	public Project(String name) {
		this.name = name;
		logger = GSLogFactory.getLog("project " + name);
		paramMap = new HashMap<ConfigKey, Object>();
	}
	
	protected Log getLogger() {
		return logger;
	}

	protected abstract ConfigKey[] getConfigKeys();

	protected ConfigKey configKeyFromName(String name) {
		for (ConfigKey key : getConfigKeys()) {
			if (key.getName().equals(name)) {
				return key;
			}
		}
		return null;
	}

	public String getName() {
		return name;
	}

	public Object configValue(ConfigKey key) {
		Object value = paramMap.get(key);
		return value == null ? key.getInfo().defaultValue() : value;
	}

	public int intConfigValue(ConfigKey key) {
		Object v = configValue(key);
		if (v instanceof Integer) {
			return (Integer) v;
		}
		throw new IllegalStateException("Value not of type int for config param " + key + ".");
	}

	public long longConfigValue(ConfigKey key) {
		Object v = configValue(key);
		if (v instanceof Long) {
			return (Long) v;
		}
		throw new IllegalStateException("Value not of type long for config param " + key + ".");
	}

	public boolean booleanConfigValue(ConfigKey key) {
		Object v = configValue(key);
		if (v instanceof Boolean) {
			return (Boolean) v;
		}
		throw new IllegalStateException("Value not of type boolean for config param " + key + ".");
	}

	public double doubleConfigValue(ConfigKey key) {
		Object v = configValue(key);
		if (v instanceof Double) {
			return (Double) v;
		}
		throw new IllegalStateException("Value not of type double for config param " + key + ".");
	}

	public String stringConfigValue(ConfigKey key) {
		Object v = configValue(key);
		if (v instanceof String) {
			return (String) v;
		}
		throw new IllegalStateException("Value not of type String for config param " + key + ".");
	}

	protected void initConfigParams(Properties... propsArray) {
		for (ConfigKey key : getConfigKeys()) {
			initConfigParam(key, propsArray);
		}
	}

	protected void initConfigParam(ConfigKey key, Properties... propsArray) {
		for (Properties properties : propsArray) {
			if (properties.getProperty(key.getName()) != null) {
				if (initConfigParam(key, key.getInfo().fromString(properties.getProperty(key.getName())))) {
					return;
				}
			}
		}
		if (!initConfigParam(key, key.getInfo().defaultValue())) {
			throw new IllegalStateException("Illegal default value " + key.getInfo().defaultValue());
		}
	}

	public boolean initConfigParam(ConfigKey key, Object value) {
		if (key.getInfo().isValueInRange(value)) {
			paramMap.put(key, value);
			return true;
		}
		return false;
	}

	protected void checkConfigProperties(Properties properties, GoalKey forGoal) {
		for (Object key : properties.keySet()) {
			ConfigKey param = configKeyFromName((String) key);
			if (param == null) {
				getLogger().warn("Unknown key " + key + ".");
			} else if (!param.isForGoal(forGoal)) {
				getLogger().warn("Bad key " + key + " for goal " + forGoal + ".");
			} else {
				String value = properties.getProperty((String) key);
				if (!param.getInfo().isValid(value)) {
					getLogger().warn("Invalid value '" + value + "' for key " + key + ".");
				} else if (!param.getInfo().isInRange(value)) {
					getLogger().warn("Value '" + value + "' out of range for key " + key + ".");
				}
			}
		}
	}
	
	protected void logParamMap() {
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Effective param map: " + paramMap);
		}
	}
}
