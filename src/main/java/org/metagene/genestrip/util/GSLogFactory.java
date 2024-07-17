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

public class GSLogFactory {
	public static Log getLog(String name) {
		return getInstance().getLogByName(name);
	}

	public static void resetN() {
		getInstance().resetNesting();
	}

	public static void incN() {
		getInstance().incNesting();
	}

	public static void decN() {
		getInstance().decNesting();
	}

	private static GSLogFactory instance;

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

	public GSLogFactory() {
		loggersByName = new HashMap<String, GSLog>();
		threadedLogOut = new ThreadLocal<PrintStream>();
		logOut = System.err;
		nestingCounter = 0;
		currentLogLevel = SimpleLog.LOG_LEVEL_INFO;
	}

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

	protected void setLogLevel(int currentLogLevel) {
		this.currentLogLevel = currentLogLevel;
	}

	public int getLogLevel() {
		return currentLogLevel;
	}

	public void resetNesting() {
		nestingCounter = 0;
	}

	public void incNesting() {
		nestingCounter++;
	}

	public void decNesting() {
		nestingCounter--;
	}

	protected void writeNesting(PrintStream logOut) {
		for (int i = 0; i < nestingCounter; i++) {
			logOut.print('>');
		}
	}

	public void setLogOut(PrintStream logOut) {
		this.logOut = logOut;
	}

	public void setLogOutForThread(PrintStream logOut) {
		threadedLogOut.set(logOut);
	}

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

	public class GSLog extends MySimpleLog {
		private static final long serialVersionUID = 1L;

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

		protected boolean isLevelEnabled(int logLevel) {
			return logLevel >= GSLogFactory.this.currentLogLevel;
		}

		@Override
		protected void write(StringBuffer buffer) {
			PrintStream tLogOut = threadedLogOut.get();
			if (tLogOut == null) {
				tLogOut = logOut;
			}
			writeNesting(tLogOut);
			tLogOut.println(buffer);
		}
	}
}
