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
package org.metagene.genestrip.store;

import org.metagene.genestrip.GSProject;

import java.io.InvalidClassException;
import java.util.Properties;

public class InvalidDatabaseClassException extends InvalidClassException {
    private static final long serialVersionUID = 1L;

    private final Properties configInfo;
    private final InvalidClassException exception;

    public InvalidDatabaseClassException(String cname, String reason, Properties configInfo, InvalidClassException e) {
        super(cname, reason);
        this.configInfo = configInfo;
        this.exception = e;
    }

    public Properties getConfigInfo() {
        return configInfo;
    }

    public InvalidClassException getException() {
        return exception;
    }

    public RuntimeException toRuntimeException() {
        String dbVersion = null;
        if (configInfo != null) {
            dbVersion = configInfo.getProperty(GSProject.GENESTRIP_VERSION);
        }
        if (dbVersion == null) {
            dbVersion = "unknown";
        }

        String title = null;
        if (configInfo != null) {
            title = configInfo.getProperty(GSProject.GENESTRIP_TITLE);
        }
        if (title == null) {
            title = "unknown";
        }

        String runtimeVersion = GSProject.getGenestripRuntimeVersion();
        if (runtimeVersion == null) {
            runtimeVersion = "unknown";
        }
        String runtimeTitle = GSProject.getGenestripRuntimeTitle();
        if (runtimeTitle == null) {
            runtimeTitle = "unknown";
        }
        return new RuntimeException("DB creation library (" + title + " version " + dbVersion + ") does not match runtime library (" +  runtimeTitle + " version " + runtimeVersion +  ").", this);
    }
}
