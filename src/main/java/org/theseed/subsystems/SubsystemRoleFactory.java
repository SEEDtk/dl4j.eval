/**
 *
 */
package org.theseed.subsystems;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.proteins.RoleMap;

/**
 * This class builds a role table from a compressed subsystem directory.  The subsystem directory is presumed to be a .tgz file.  We
 * don't decompress it because subsystems are not designed to be machine-portable.  Instead, we process the files in place to determine
 * which roles are in valid subsystems.
 *
 * A GZipCompressorInputStream converts the compressed file into an archive.  The TarArchiveInputStream then interprets the archive and
 * returns directories and files.  Each directory is a subsystem.  Inside the subsystem, we care about the CLASSIFICATION and
 * EXCHANGEABLE files (which determine whether or not the spreadsheet is valid) and the "spreadsheet" file, which contains the actual
 * roles.  Any roles found will be added to the role table.
 *
 * @author Bruce Parrello
 *
 */
public class SubsystemRoleFactory {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(SubsystemRoleFactory.class);
    /** role table being built */
    private RoleMap roleTable;
    /** list of roles in the current subsystem */
    List<String> currRoles;
    /** TRUE if this subsystem is public (good) */
    private boolean currPublic;
    /** TRUE if this subsystem is experimental (bad) */
    private boolean currExperimental;
    /** name of the current subsystem */
    private String currName;

    /**
     * Initialize a subsystem role factory.
     */
    public SubsystemRoleFactory() {
        this.roleTable = new RoleMap();
    }

    /**
     * Initialize a subsystem role factory, pre-loading the specified roles.
     *
     * @param initRoles		role-map save file containing the roles to pre-load
     */
    public SubsystemRoleFactory(File initRoles) {
        this.roleTable = RoleMap.load(initRoles);
    }

    /**
     * Read a subsystem archive and create a table of the subsystem roles.
     *
     * @param subsystemArchive	name of the subsystem archive file
     *
     * @return the table of subsystem roles
     *
     * @throws IOException
     */
    public static RoleMap processArchive(File subsystemArchive) throws IOException {
        SubsystemRoleFactory factory = new SubsystemRoleFactory();
        factory.fillSubsystems(subsystemArchive);
        // Return the role table.
        return factory.roleTable;
    }

    /**
     * Read a subsystem archive and create a table of the subsystem roles.
     *
     * @param subsystemArchive	name of the subsystem archive file
     * @param initRoles			saved roles to use to initialize the role table
     *
     * @return the table of subsystem roles
     *
     * @throws IOException
     */
    public static RoleMap processArchive(File subsystemArchive, File initRoles) throws IOException {
        SubsystemRoleFactory factory = new SubsystemRoleFactory(initRoles);
        factory.fillSubsystems(subsystemArchive);
        // Return the role table.
        return factory.roleTable;
    }

    /**
     * Fill this subsystem role factory with the roles from the specified subsystem archive.
     *
     * @param subsystemArchive	subsystem archive file to read
     *
     * @throws IOException
     */
    protected void fillSubsystems(File subsystemArchive) throws IOException {
        // Open up the input stream.
        try (FileInputStream fileStream = new FileInputStream(subsystemArchive)) {
            GzipCompressorInputStream gzIn = new GzipCompressorInputStream(fileStream);
            TarArchiveInputStream tarIn = new TarArchiveInputStream(gzIn);
            // Allocate the role string list.
            this.currRoles = new ArrayList<String>(50);
            // Now we loop through the archive entries.  A key point here is that every entry name is fully-qualified
            // and begins with "Subsystems/".  We use this information to insure we don't recurse into subdirectories.
            TarArchiveEntry curr = tarIn.getNextTarEntry();
            while (curr != null) {
                if (curr.isDirectory()) {
                    // Insure this is really a subsystem.
                    String name = StringUtils.substring(curr.getName(), 11, -1);
                    if (! StringUtils.contains(name, "/")) {
                        // Here we have a new subsystem.  Check the old one.
                        this.closeSubsystem();
                        // Initialize the new one.
                        this.openSubsystem(name);
                    }
                } else {
                    // Here we have a subsystem file.  Check the name.
                    String name = StringUtils.substring(curr.getName(), 12 + this.currName.length());
                    switch (name) {
                    case "CLASSIFICATION" :
                        this.checkExperimental(tarIn);
                        break;
                    case "EXCHANGEABLE" :
                    case "EXCHANGABLE"  :
                        // Here the subsystem is public.  Set the public bit.
                        this.currPublic = true;
                        break;
                    case "spreadsheet" :
                        // Here we must read the roles.
                        this.readRoles(tarIn);
                    }
                }
                // Get the next archive entry.
                curr = tarIn.getNextTarEntry();
            }
        }
        // Close the last subsystem.
        this.closeSubsystem();
    }



    /**
     * Check the classification data for an experimental subsystem.
     *
     * @param tarIn		current input stream
     *
     * @throws IOException
     */
    private void checkExperimental(TarArchiveInputStream tarIn) throws IOException {
        List<String> lines = IOUtils.readLines(tarIn, Charset.defaultCharset());
        if (lines.size() > 0 && StringUtils.containsIgnoreCase(lines.get(0), "experimental"))
            this.currExperimental = true;
    }



    /**
     * This method reads roles from the current file.
     *
     * @param tarIn		current archive stream for the subsystem directory
     *
     * @throws IOException
     */
    private void readRoles(TarArchiveInputStream tarIn) throws IOException {
        List<String> lines = IOUtils.readLines(tarIn, Charset.defaultCharset());
        if (lines.size() > 0) {
            Iterator<String> iter = lines.iterator();
            String line = iter.next();
            while (! line.contentEquals("//")) {
                String roleInfo = StringUtils.substringAfter(line, "\t");
                String[] parts = Feature.rolesOfFunction(roleInfo);
                for (String part : parts)
                    this.currRoles.add(part);
                if (iter.hasNext())
                    line = iter.next();
                else
                    line = "//";
            }
        }
    }

    /**
     * Initialize for processing a new subsystem.
     *
     * @param name		name of this subsystem.
     */
    private void openSubsystem(String name) {
        // Denote no roles are queued.
        this.currRoles.clear();
        // Denote the subsystem is not public or experimental.
        this.currExperimental = false;
        this.currPublic = false;
        // Save the name.
        this.currName = name;
    }

    /**
     * Process this subsystem.
     */
    private void closeSubsystem() {
        // Only proceed if we have a name.
        if (this.currName != null) {
            if (! this.currPublic)
                log.info("Private subsystem {} skipped.", this.currName);
            else if (this.currExperimental)
                log.info("Experimental subsystem {} skipped.", this.currName);
            else {
                // Here the subsystem is valid.
                log.info("{} roles found in subsystem {}.", this.currRoles.size(), this.currName);
                for (String role : this.currRoles)
                    this.roleTable.findOrInsert(role);
            }
        }
    }

}
