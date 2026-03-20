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
