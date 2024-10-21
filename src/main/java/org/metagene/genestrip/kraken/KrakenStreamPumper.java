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
package org.metagene.genestrip.kraken;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import org.apache.commons.exec.util.DebugUtils;

/**
 * Copies all data from an input stream to an output stream.
 *
 * @version $Id: StreamPumper.java 1557263 2014-01-10 21:18:09Z ggregory $
 */
public class KrakenStreamPumper implements Runnable {

	/** the default size of the internal buffer for copying the streams */
	private static final int DEFAULT_SIZE = 1024;

	/** the input stream to pump from */
	private final InputStream is;

	/** the output stream to pmp into */
	private final OutputStream os;

	/** the size of the internal buffer for copying the streams */
	private final int size;

	/** was the end of the stream reached */
	private boolean finished;

	/** close the output stream when exhausted */
	private final boolean closeWhenExhausted;

	/**
	 * Create a new stream pumper.
	 * 
	 * @param is                 input stream to read data from
	 * @param os                 output stream to write data to.
	 * @param closeWhenExhausted if true, the output stream will be closed when the
	 *                           input is exhausted.
	 */
	public KrakenStreamPumper(final InputStream is, final OutputStream os, final boolean closeWhenExhausted) {
		this.is = is;
		this.os = os;
		this.size = DEFAULT_SIZE;
		this.closeWhenExhausted = closeWhenExhausted;
	}

	/**
	 * Create a new stream pumper.
	 *
	 * @param is                 input stream to read data from
	 * @param os                 output stream to write data to.
	 * @param closeWhenExhausted if true, the output stream will be closed when the
	 *                           input is exhausted.
	 * @param size               the size of the internal buffer for copying the
	 *                           streams
	 */
	public KrakenStreamPumper(final InputStream is, final OutputStream os, final boolean closeWhenExhausted,
			final int size) {
		this.is = is;
		this.os = os;
		this.size = size > 0 ? size : DEFAULT_SIZE;
		this.closeWhenExhausted = closeWhenExhausted;
	}

	/**
	 * Create a new stream pumper.
	 * 
	 * @param is input stream to read data from
	 * @param os output stream to write data to.
	 */
	public KrakenStreamPumper(final InputStream is, final OutputStream os) {
		this(is, os, false);
	}

	/**
	 * Copies data from the input stream to the output stream. Terminates as soon as
	 * the input stream is closed or an error occurs.
	 */
	public void run() {
		synchronized (this) {
			// Just in case this object is reused in the future
			finished = false;
		}

		try {
			doWork(is, os);
		} catch (final Exception e) {
			// nothing to do - happens quite often with watchdog
		} finally {
			if (closeWhenExhausted && os != null) {
				try {
					os.close();
				} catch (final IOException e) {
					final String msg = "Got exception while closing exhausted output stream";
					DebugUtils.handleException(msg, e);
				}
			}
			synchronized (this) {
				finished = true;
				notifyAll();
			}
		}
	}

	protected void doWork(InputStream is, OutputStream os) throws IOException {
		if (os != null) {
			int length;
			final byte[] buf = new byte[this.size];
			while ((length = is.read(buf)) > 0) {
				os.write(buf, 0, length);
			}
		}
	}

	/**
	 * Tells whether the end of the stream has been reached.
	 * 
	 * @return true is the stream has been exhausted.
	 */
	public synchronized boolean isFinished() {
		return finished;
	}

	/**
	 * This method blocks until the stream pumper finishes.
	 * 
	 * @exception InterruptedException if any thread interrupted the current thread
	 *                                 before or while the current thread was
	 *                                 waiting for a notification.
	 * @see #isFinished()
	 */
	public synchronized void waitFor() throws InterruptedException {
		while (!isFinished()) {
			wait();
		}
	}
}
