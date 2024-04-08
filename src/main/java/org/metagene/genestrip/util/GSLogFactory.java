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

	public GSLogFactory() {
		loggersByName = new HashMap<String, GSLog>();
		threadedLogOut = new ThreadLocal<PrintStream>();
		logOut = System.err;
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

	public class GSLog extends SimpleLog {
		private static final long serialVersionUID = 1L;

		public GSLog(String name) {
			super(name);
		}

		@Override
		protected void write(StringBuffer buffer) {
			PrintStream tLogOut = threadedLogOut.get();
			if (tLogOut != null) {
				tLogOut.println(buffer);
			} else {
				logOut.println(buffer);
			}
		}
	}
}
