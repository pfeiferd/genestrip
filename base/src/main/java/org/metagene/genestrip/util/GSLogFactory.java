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
package org.metagene.genestrip.util;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.commons.logging.impl.SimpleLog;

/**
 * A logging facility built on Apache Commons Logging's {@link SimpleLog} that shares a single,
 * dynamically adjustable log level across all loggers, supports redirecting output per thread, and
 * prefixes each message with {@code '>'} characters to reflect the current make nesting depth. A
 * lazily created singleton backs the static convenience methods.
 */
public class GSLogFactory {
	/**
	 * Returns the shared logger with the given name, creating it if necessary.
	 *
	 * @param name the logger name
	 * @return the shared logger with the given name
	 */
	public static Log getLog(String name) {
		return getInstance().getLogByName(name);
	}

	/**
	 * Resets the make-nesting indentation level of the shared factory to zero.
	 */
	public static void resetN() {
		getInstance().resetNesting();
	}

	/**
	 * Increases the make-nesting indentation level of the shared factory by one.
	 */
	public static void incN() {
		getInstance().incNesting();
	}

	/**
	 * Decreases the make-nesting indentation level of the shared factory by one.
	 */
	public static void decN() {
		getInstance().decNesting();
	}

	private static GSLogFactory instance;

	/**
	 * Returns the lazily created singleton factory.
	 *
	 * @return the singleton factory instance
	 */
	public static GSLogFactory getInstance() {
		if (instance == null) {
			instance = new GSLogFactory();
		}
		return instance;
	}

	private final Map<String, GSLog> loggersByName;
	private PrintStream logOut;
	private final ThreadLocal<PrintStream> threadedLogOut;
	private int nestingCounter;
	private volatile int currentLogLevel;

	/**
	 * Creates a factory with an empty logger cache, output to {@code System.err}, zero nesting and
	 * the default {@code INFO} log level.
	 */
	public GSLogFactory() {
		loggersByName = new HashMap<String, GSLog>();
		threadedLogOut = new ThreadLocal<PrintStream>();
		logOut = System.err;
		nestingCounter = 0;
		currentLogLevel = SimpleLog.LOG_LEVEL_INFO;
	}

	/**
	 * Sets the shared log level from its name: one of {@code all}, {@code trace}, {@code debug},
	 * {@code info}, {@code warn}, {@code error}, {@code fatal} or {@code off} (case-insensitive).
	 *
	 * @param lvl the log level name
	 * @throws IllegalArgumentException if the name is not a recognized level
	 */
	public void setLogLevel(String lvl) {
		if ("all".equalsIgnoreCase(lvl)) {
			setLogLevel(SimpleLog.LOG_LEVEL_ALL);
		} else if ("trace".equalsIgnoreCase(lvl)) {
			setLogLevel(SimpleLog.LOG_LEVEL_TRACE);
		} else if ("debug".equalsIgnoreCase(lvl)) {
			setLogLevel(SimpleLog.LOG_LEVEL_DEBUG);
		} else if ("info".equalsIgnoreCase(lvl)) {
			setLogLevel(SimpleLog.LOG_LEVEL_INFO);
		} else if ("warn".equalsIgnoreCase(lvl)) {
			setLogLevel(SimpleLog.LOG_LEVEL_WARN);
		} else if ("error".equalsIgnoreCase(lvl)) {
			setLogLevel(SimpleLog.LOG_LEVEL_ERROR);
		} else if ("fatal".equalsIgnoreCase(lvl)) {
			setLogLevel(SimpleLog.LOG_LEVEL_FATAL);
		} else if ("off".equalsIgnoreCase(lvl)) {
			setLogLevel(SimpleLog.LOG_LEVEL_OFF);
		} else {
			throw new IllegalArgumentException("Illegal log level: " + lvl);
		}		
	}

	/**
	 * Sets the shared log level to the given {@link SimpleLog} level constant.
	 *
	 * @param currentLogLevel the numeric log level to set
	 */
	protected void setLogLevel(int currentLogLevel) {
		this.currentLogLevel = currentLogLevel;
	}

	/**
	 * Returns the current shared log level.
	 *
	 * @return the current numeric log level
	 */
	public int getLogLevel() {
		return currentLogLevel;
	}

	/**
	 * Resets the nesting indentation level (the number of {@code '>'} prefixes on log lines) to zero.
	 */
	public void resetNesting() {
		nestingCounter = 0;
	}

	/**
	 * Increases the nesting indentation level (the {@code '>'} prefix depth on log lines) by one.
	 */
	public void incNesting() {
		nestingCounter++;
	}

	/**
	 * Decreases the nesting indentation level (the {@code '>'} prefix depth on log lines) by one.
	 */
	public void decNesting() {
		nestingCounter--;
	}

	/**
	 * Writes the current nesting indentation ({@code '>'} prefixes) to the given output stream.
	 *
	 * @param logOut the stream to write the indentation to
	 */
	protected void writeNesting(PrintStream logOut) {
		for (int i = 0; i < nestingCounter; i++) {
			logOut.print('>');
		}
	}

	/**
	 * Sets the default output stream to which log messages are written.
	 *
	 * @param logOut the default log output stream
	 */
	public void setLogOut(PrintStream logOut) {
		this.logOut = logOut;
	}

	/**
	 * Sets a log output stream used only by the current thread, overriding the default output stream
	 * for messages logged from that thread.
	 *
	 * @param logOut the per-thread log output stream
	 */
	public void setLogOutForThread(PrintStream logOut) {
		threadedLogOut.set(logOut);
	}

	/**
	 * Returns the shared logger with the given name, creating and caching it if necessary.
	 *
	 * @param name the logger name
	 * @return the shared logger with the given name
	 */
	public GSLog getLogByName(String name) {
		if (logOut == null) {
			LogFactory.getLog(name);
		}

		GSLog log = loggersByName.get(name);
		if (log == null) {
			log = new GSLog(name);
			loggersByName.put(name, log);
		}
		return log;
	}

	// Required to just to run the static initializer from below.
	private static class MySimpleLog extends SimpleLog {
		private static final long serialVersionUID = 1L;

		static {
			SimpleLog.showDateTime = false;
		}

		public MySimpleLog(String name) {
			super(name);
		}
	}

	/**
	 * A {@link SimpleLog} whose effective level is governed centrally by the enclosing factory and
	 * whose output honors the factory's nesting indentation and per-thread output redirection.
	 */
	public class GSLog extends MySimpleLog {
		private static final long serialVersionUID = 1L;

		/**
		 * Creates a logger with the given name backed by the enclosing factory.
		 *
		 * @param name the logger name
		 */
		public GSLog(String name) {
			super(name);
		}

		// Ensure central dynamic setting of log level.
		@Override
		public int getLevel() {
			return GSLogFactory.this.currentLogLevel;
		}

		@Override
		public void setLevel(int currentLogLevel) {
			// Do nothing on purpose.
		}

		/**
		 * Returns this logger's name.
		 *
		 * @return this logger's name
		 */
		public String getName() {
			return logName;
		}

		protected boolean isLevelEnabled(int logLevel) {
			return logLevel >= GSLogFactory.this.currentLogLevel;
		}

		@Override
		protected void write(StringBuffer buffer) {
			PrintStream tLogOut = threadedLogOut.get();
			if (tLogOut == null) {
				tLogOut = logOut;
			}
			synchronized (tLogOut) {
				writeNesting(tLogOut);
				tLogOut.println(buffer);
			}
		}
	}
}
