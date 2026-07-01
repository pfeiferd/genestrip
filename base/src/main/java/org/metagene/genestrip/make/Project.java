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
import org.metagene.genestrip.util.GSLogFactory;

/**
 * Base class for a project: a named holder of typed configuration parameters (each a {@link ConfigKey}
 * whose value is validated against its {@link ConfigParamInfo}) plus free-form additional
 * properties. Goals read their configuration through the project. Subclasses declare the supported
 * configuration keys via {@link #getConfigKeys()}.
 */
public abstract class Project {
	/** Additional-property key under which the project name is stored. */
	public static final String PROJECT_NAME = "projectName";

	private final Log logger;
	private final String name;
	private final Map<ConfigKey, Object> paramMap;
	private final Properties additionalProperties;
	
	/**
	 * Creates a project with the given name, which is also recorded as the {@value #PROJECT_NAME}
	 * additional property.
	 *
	 * @param name the project name
	 */
	public Project(String name) {
		this.name = name;
		logger = GSLogFactory.getLog("project " + name);
		paramMap = new HashMap<ConfigKey, Object>();
		additionalProperties = new Properties();
		setAdditionalProperty(PROJECT_NAME, name);
	}

	/**
	 * Sets a free-form additional property.
	 *
	 * @param key the property key
	 * @param value the property value
	 */
	public void setAdditionalProperty(String key, String value) {
		additionalProperties.setProperty(key, String.valueOf(value));
	}

	/**
	 * Returns the value of a free-form additional property.
	 *
	 * @param key the property key
	 * @return the property value, or {@code null} if unset
	 */
	public String getAdditionalProperty(String key) {
		return additionalProperties.getProperty(key);
	}

	/**
	 * Returns a copy of the free-form additional properties.
	 *
	 * @return a copy of the additional properties
	 */
	public Properties getAdditionalProperties() {
		Properties props = new Properties();
		props.putAll(additionalProperties);
		return props;
	}

	/**
	 * Returns the configured parameter values as properties keyed by config key name.
	 *
	 * @return the config parameters as properties
	 */
	public Properties getConfigAsProperties() {
		Properties props = new Properties();
		for (ConfigKey key : paramMap.keySet()) {
			props.put(key.getName(), String.valueOf(paramMap.get(key)));
		}
		return props;
	}

	/**
	 * Returns the config parameters and the additional properties combined into a single
	 * {@link Properties}.
	 *
	 * @return the config parameters and additional properties combined
	 */
	public Properties getAllAsProperties() {
		Properties props = new Properties();
		props.putAll(getConfigAsProperties());
		props.putAll(additionalProperties);
		return props;
	}
	
	/**
	 * Returns the logger of this project.
	 *
	 * @return the project logger
	 */
	protected Log getLogger() {
		return logger;
	}

	/**
	 * Returns all configuration keys supported by this project.
	 *
	 * @return the supported configuration keys
	 */
	protected abstract ConfigKey[] getConfigKeys();

	/**
	 * Returns the supported config key with the given name, or {@code null} if there is none.
	 *
	 * @param name the config key name to look up
	 * @return the matching config key, or {@code null} if there is none
	 */
	protected ConfigKey configKeyFromName(String name) {
		for (ConfigKey key : getConfigKeys()) {
			if (key.getName().equals(name)) {
				return key;
			}
		}
		return null;
	}

	/**
	 * Returns the project name.
	 *
	 * @return the project name
	 */
	public String getName() {
		return name;
	}

	/**
	 * Returns the configured value for the given key, or the key's default value if it is unset.
	 *
	 * @param key the config key
	 * @return the configured value, or the key's default value if unset
	 */
	public Object configValue(ConfigKey key) {
		Object value = paramMap.get(key);
		return value == null ? key.getInfo().defaultValue() : value;
	}

	/**
	 * Returns the value of the given key as an {@code int}.
	 *
	 * @param key the config key
	 * @return the value as an {@code int}
	 * @throws IllegalStateException if the value is not an {@code int}
	 */
	public int intConfigValue(ConfigKey key) {
		Object v = configValue(key);
		if (v instanceof Integer) {
			return (Integer) v;
		}
		throw new IllegalStateException("Value not of type int for config param " + key + ".");
	}

	/**
	 * Returns the value of the given key as a {@code long}.
	 *
	 * @param key the config key
	 * @return the value as a {@code long}
	 * @throws IllegalStateException if the value is not a {@code long}
	 */
	public long longConfigValue(ConfigKey key) {
		Object v = configValue(key);
		if (v instanceof Long) {
			return (Long) v;
		}
		throw new IllegalStateException("Value not of type long for config param " + key + ".");
	}

	/**
	 * Returns the value of the given key as a {@code boolean}.
	 *
	 * @param key the config key
	 * @return the value as a {@code boolean}
	 * @throws IllegalStateException if the value is not a {@code boolean}
	 */
	public boolean booleanConfigValue(ConfigKey key) {
		Object v = configValue(key);
		if (v instanceof Boolean) {
			return (Boolean) v;
		}
		throw new IllegalStateException("Value not of type boolean for config param " + key + ".");
	}

	/**
	 * Returns the value of the given key as a {@code double}.
	 *
	 * @param key the config key
	 * @return the value as a {@code double}
	 * @throws IllegalStateException if the value is not a {@code double}
	 */
	public double doubleConfigValue(ConfigKey key) {
		Object v = configValue(key);
		if (v instanceof Double) {
			return (Double) v;
		}
		throw new IllegalStateException("Value not of type double for config param " + key + ".");
	}

	/**
	 * Returns the value of the given key as a {@link String}.
	 *
	 * @param key the config key
	 * @return the value as a {@link String}
	 * @throws IllegalStateException if the value is not a {@link String}
	 */
	public String stringConfigValue(ConfigKey key) {
		Object v = configValue(key);
		if (v instanceof String) {
			return (String) v;
		}
		throw new IllegalStateException("Value not of type String for config param " + key + ".");
	}

	/**
	 * Initializes all supported config parameters from the given property sources.
	 *
	 * @param propsArray the property sources, in decreasing priority
	 */
	protected void initConfigParams(Properties... propsArray) {
		for (ConfigKey key : getConfigKeys()) {
			initConfigParam(key, propsArray);
		}
	}

	/**
	 * Initializes a single config parameter from the first property source that supplies a value the
	 * key accepts, falling back to the key's default value.
	 *
	 * @param key the config key to initialize
	 * @param propsArray the property sources, in decreasing priority
	 * @throws IllegalStateException if even the default value is illegal
	 */
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

	/**
	 * Sets the value of the given config parameter if it is within the key's allowed range.
	 *
	 * @param key the config key to set
	 * @param value the value to store
	 * @return whether the value was accepted and stored
	 */
	public boolean initConfigParam(ConfigKey key, Object value) {
		if (key.getInfo().isValueInRange(value)) {
			paramMap.put(key, value);
			return true;
		}
		return false;
	}

	/**
	 * Logs warnings for entries of the given properties that reference an unknown key, a key that is
	 * not valid for the given goal, or an invalid or out-of-range value.
	 *
	 * @param properties the properties to check
	 * @param forGoal the goal the properties are checked against
	 */
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
	
	/**
	 * Logs the effective configuration parameter map at info level.
	 */
	public void logParamMap() {
		if (getLogger().isInfoEnabled()) {
			getLogger().info("Effective param map: " + paramMap);
		}
	}
}
