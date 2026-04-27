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

public abstract class ConfigParamInfo<V> {
	private V defaultValue;

	public ConfigParamInfo(V defaultValue) {
		this.defaultValue = defaultValue;
	}

	public V defaultValue() {
		return defaultValue;
	}

	public String getMDDefaultValue() {
		return defaultValue == null ? "" : defaultValue.toString();
	}

	public boolean isValid(String s) {
		return fromString(s) != null;
	}

	public boolean isInRange(String s) {
		V v = fromString(s);
		return v != null && isValueInRange(v);
	}

	public boolean isValueInRange(Object o) {
		return o != null;
	}

	public String getMDRangeDescriptor() {
		return "";
	}

	public String getTypeDescriptor() {
		return "";
	}

	protected abstract V fromString(String s);

	public static class IntConfigParamInfo extends ConfigParamInfo<Integer> {
		private final int min;
		private final int max;

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

	public static class LongConfigParamInfo extends ConfigParamInfo<Long> {
		private final long min;
		private final long max;

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

	public static class DoubleConfigParamInfo extends ConfigParamInfo<Double> {
		private final double min;
		private final double max;

		public DoubleConfigParamInfo(double min, double max, double defaultValue) {
			super(defaultValue);
			this.min = min;
			this.max = max;
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
				Double d = (Double) value;
				return d >= min && d <= max;
			}
			return false;
		}

		@Override
		public String getMDRangeDescriptor() {
			return "[" + min + ", " + max + "]";
		}

		@Override
		public String getTypeDescriptor() {
			return "double";
		}
	}

	public static class BooleanConfigParamInfo extends ConfigParamInfo<Boolean> {
		public BooleanConfigParamInfo(boolean defaultValue) {
			super(defaultValue);
		}

		@Override
		protected Boolean fromString(String s) {
			return s == null ? null : Boolean.valueOf(s.trim());
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

	public static class StringConfigParamInfo extends ConfigParamInfo<String> {
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

	public static abstract class ListConfigParamInfo<E> extends ConfigParamInfo<List<E>> {
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
