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

import java.util.List;

/**
 * Describes a configuration parameter's value type: its default value, how to parse a value from a
 * string and how to check that a value is within the allowed range. Concrete subclasses exist for
 * {@code int}, {@code long}, {@code double}, {@code boolean}, {@link String} and list parameters.
 * The {@code getMD...} methods provide descriptors used when generating Markdown documentation.
 *
 * @param <V> the value type of the parameter
 */
public abstract class ConfigParamInfo<V> {
	private V defaultValue;

	/**
	 * Creates a descriptor with the given default value.
	 *
	 * @param defaultValue the default value of the parameter
	 */
	public ConfigParamInfo(V defaultValue) {
		this.defaultValue = defaultValue;
	}

	/**
	 * Returns the default value of this parameter.
	 *
	 * @return the default value
	 */
	public V defaultValue() {
		return defaultValue;
	}

	/**
	 * Returns the default value rendered for Markdown documentation.
	 *
	 * @return the default value as a Markdown string
	 */
	public String getMDDefaultValue() {
		return defaultValue == null ? "" : defaultValue.toString();
	}

	/**
	 * Returns whether the given string can be parsed into a value of this parameter's type.
	 *
	 * @param s the string to parse
	 * @return {@code true} if the string parses to a valid value
	 */
	public boolean isValid(String s) {
		return fromString(s) != null;
	}

	/**
	 * Returns whether the given string parses to a value that is within the allowed range.
	 *
	 * @param s the string to parse
	 * @return {@code true} if the parsed value is within the allowed range
	 */
	public boolean isInRange(String s) {
		V v = fromString(s);
		return v != null && isValueInRange(v);
	}

	/**
	 * Returns whether the given value is within the allowed range; by default any non-{@code null}
	 * value.
	 *
	 * @param o the value to check
	 * @return {@code true} if the value is within the allowed range
	 */
	public boolean isValueInRange(Object o) {
		return o != null;
	}

	/**
	 * Returns a Markdown description of the allowed range; empty by default.
	 *
	 * @return the range description as a Markdown string
	 */
	public String getMDRangeDescriptor() {
		return "";
	}

	/**
	 * Returns a short name of the value type (e.g. {@code "int"}); empty by default.
	 *
	 * @return the value type descriptor
	 */
	public String getTypeDescriptor() {
		return "";
	}

	/**
	 * Parses a value from the given string, returning {@code null} if it cannot be parsed.
	 *
	 * @param s the string to parse
	 * @return the parsed value, or {@code null} if it cannot be parsed
	 */
	protected abstract V fromString(String s);

	/**
	 * A {@link ConfigParamInfo} for {@code int} parameters bounded by an inclusive {@code [min, max]}
	 * range.
	 */
	public static class IntConfigParamInfo extends ConfigParamInfo<Integer> {
		private final int min;
		private final int max;

		/**
		 * Creates an {@code int} parameter descriptor.
		 *
		 * @param min          the inclusive lower bound
		 * @param max          the inclusive upper bound
		 * @param defaultValue the default value
		 */
		public IntConfigParamInfo(int min, int max, int defaultValue) {
			super(defaultValue);
			this.min = min;
			this.max = max;
		}

		@Override
		protected Integer fromString(String s) {
			try {
				return s == null ? null : Integer.valueOf(s.trim());
			} catch (NumberFormatException e) {
				return null;
			}
		}

		@Override
		public boolean isValueInRange(Object value) {
			if (value instanceof Integer) {
				Integer i = (Integer) value;
				return i >= min && i <= max;
			}
			return false;
		}

		@Override
		public String getMDRangeDescriptor() {
			return "[" + min + ", " + max + "]";
		}

		@Override
		public String getTypeDescriptor() {
			return "int";
		}
	}

	/**
	 * A {@link ConfigParamInfo} for {@code long} parameters bounded by an inclusive {@code [min, max]}
	 * range.
	 */
	public static class LongConfigParamInfo extends ConfigParamInfo<Long> {
		private final long min;
		private final long max;

		/**
		 * Creates a {@code long} parameter descriptor.
		 *
		 * @param min          the inclusive lower bound
		 * @param max          the inclusive upper bound
		 * @param defaultValue the default value
		 */
		public LongConfigParamInfo(long min, long max, long defaultValue) {
			super(defaultValue);
			this.min = min;
			this.max = max;
		}

		@Override
		protected Long fromString(String s) {
			try {
				return s == null ? null : Long.valueOf(s.trim());
			} catch (NumberFormatException e) {
				return null;
			}
		}

		@Override
		public boolean isValueInRange(Object value) {
			if (value instanceof Long) {
				Long i = (Long) value;
				return i >= min && i <= max;
			}
			return false;
		}

		@Override
		public String getMDRangeDescriptor() {
			return "[" + min + ", " + max + "]";
		}

		@Override
		public String getTypeDescriptor() {
			return "long";
		}
	}

	/**
	 * A {@link ConfigParamInfo} for {@code double} parameters bounded by the closed interval
	 * {@code [min, max]}, or the open interval {@code (min, max)} when created as exclusive.
	 */
	public static class DoubleConfigParamInfo extends ConfigParamInfo<Double> {
		private final double min;
		private final double max;
		// When true the value must lie in the OPEN interval (min, max) rather than [min, max].
		private final boolean exclusive;

		/**
		 * Creates a {@code double} parameter descriptor with a closed range.
		 *
		 * @param min          the inclusive lower bound
		 * @param max          the inclusive upper bound
		 * @param defaultValue the default value
		 */
		public DoubleConfigParamInfo(double min, double max, double defaultValue) {
			this(min, max, defaultValue, false);
		}

		/**
		 * Creates a {@code double} parameter descriptor with an optionally exclusive range.
		 *
		 * @param min          the lower bound
		 * @param max          the upper bound
		 * @param defaultValue the default value
		 * @param exclusive    {@code true} for an open interval {@code (min, max)}, {@code false} for a
		 *                     closed interval {@code [min, max]}
		 */
		public DoubleConfigParamInfo(double min, double max, double defaultValue, boolean exclusive) {
			super(defaultValue);
			this.min = min;
			this.max = max;
			this.exclusive = exclusive;
		}

		@Override
		protected Double fromString(String s) {
			try {
				return s == null ? null : Double.valueOf(s.trim());
			} catch (NumberFormatException e) {
				return null;
			}
		}

		@Override
		public boolean isValueInRange(Object value) {
			if (value instanceof Double) {
				double d = (Double) value;
				return exclusive ? (d > min && d < max) : (d >= min && d <= max);
			}
			return false;
		}

		@Override
		public String getMDRangeDescriptor() {
			return (exclusive ? "(" : "[") + min + ", " + max + (exclusive ? ")" : "]");
		}

		@Override
		public String getTypeDescriptor() {
			return "double";
		}
	}

	/**
	 * A {@link ConfigParamInfo} for {@code boolean} parameters, accepting only {@code "true"} and
	 * {@code "false"} (case-insensitive).
	 */
	public static class BooleanConfigParamInfo extends ConfigParamInfo<Boolean> {
		/**
		 * Creates a {@code boolean} parameter descriptor.
		 *
		 * @param defaultValue the default value
		 */
		public BooleanConfigParamInfo(boolean defaultValue) {
			super(defaultValue);
		}

		@Override
		protected Boolean fromString(String s) {
			if (s == null) {
				return null;
			}
			String t = s.trim();
			// Return null (not false) for anything that is not "true"/"false" so that a typo is
			// reported as an invalid value instead of being silently coerced to false.
			if (t.equalsIgnoreCase("true")) {
				return Boolean.TRUE;
			}
			if (t.equalsIgnoreCase("false")) {
				return Boolean.FALSE;
			}
			return null;
		}

		@Override
		public String getMDRangeDescriptor() {
			return "";
		}

		@Override
		public String getTypeDescriptor() {
			return "boolean";
		}
	}

	/**
	 * A {@link ConfigParamInfo} for {@link String} parameters; values are trimmed when parsed.
	 */
	public static class StringConfigParamInfo extends ConfigParamInfo<String> {
		/**
		 * Creates a {@link String} parameter descriptor.
		 *
		 * @param defaultValue the default value
		 */
		public StringConfigParamInfo(String defaultValue) {
			super(defaultValue);
		}

		@Override
		protected String fromString(String s) {
			return s == null ? null : s.trim();
		}

		@Override
		public String getMDRangeDescriptor() {
			return "";
		}

		@Override
		public String getTypeDescriptor() {
			return "String";
		}
	}

	/**
	 * A {@link ConfigParamInfo} for list-valued parameters whose elements are of type {@code E}; the
	 * Markdown default value is rendered as the comma-separated elements.
	 *
	 * @param <E> the element type of the list
	 */
	public static abstract class ListConfigParamInfo<E> extends ConfigParamInfo<List<E>> {
		/**
		 * Creates a list parameter descriptor.
		 *
		 * @param defaultValue the default list value
		 */
		public ListConfigParamInfo(List<E> defaultValue) {
			super(defaultValue);
		}

		public String getMDDefaultValue() {
			if (defaultValue() == null) {
				return "";
			}
			StringBuilder sb = new StringBuilder();
			boolean first = true;
			for (E e : defaultValue()) {
				if (!first) {
					sb.append(',');
				}
				first = false;
				sb.append(e);
			}
			return sb.toString();
		}
	}
}
