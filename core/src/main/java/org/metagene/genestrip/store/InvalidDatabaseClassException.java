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

    public static RuntimeException convertToRuntimeException(InvalidClassException e) {
        if (e instanceof InvalidDatabaseClassException) {
            InvalidDatabaseClassException e2 = (InvalidDatabaseClassException) e;
            Properties properties = e2.getConfigInfo();
            String dbVersion = properties.getProperty(GSProject.GENESTRIP_DB_VERSION);
            if (dbVersion == null) {
                dbVersion = "unknown";
            }
            String runtimeVersion = GSProject.getGenestripRuntimeVersion();
            return new RuntimeException("Database file version (" + dbVersion + ") does not match Genestrip library version (" + runtimeVersion + ").", e);
        }
        else {
            return new RuntimeException("Database file version does not match Genestrip library version.", e);
        }
    }

    private final Properties configInfo;
    private final String className;

    public InvalidDatabaseClassException(String cname, String reason, Properties configInfo) {
        super(null, reason);
        this.className = cname;
        this.configInfo = configInfo;
    }

    public String getClassName() {
        return className;
    }

    public Properties getConfigInfo() {
        return configInfo;
    }
}
