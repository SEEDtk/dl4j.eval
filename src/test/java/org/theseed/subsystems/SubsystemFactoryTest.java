/**
 *
 */
package org.theseed.subsystems;

import org.junit.jupiter.api.Test;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.Collection;

import org.theseed.io.LineReader;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;

/**
 * @author Bruce Parrello
 *
 */
public class SubsystemFactoryTest {

    @Test
    public void testFactory() throws IOException {
        File subFile = new File("data", "subsystems.tgz");
        RoleMap subRoles = SubsystemRoleFactory.processArchive(subFile);
        File testFile = new File("data", "roles.in.subs.txt");
        Collection<String> roles = LineReader.readList(testFile);
        assertThat(subRoles.size(), equalTo(roles.size()));
        for (String roleDesc : roles) {
            Role role = subRoles.getByName(roleDesc);
            assertThat(roleDesc, role, not(nullValue()));
        }
    }

    @Test
    public void testPreLoad() throws IOException {
        File subFile = new File("data", "subsystems.tgz");
        File initFile = new File("data", "roles.init.tbl");
        RoleMap subRoles = SubsystemRoleFactory.processArchive(subFile, initFile);
        File testFile = new File("data", "roles.in.subs.txt");
        Collection<String> roles = LineReader.readList(testFile);
        assertThat(subRoles.size(), equalTo(roles.size() + 31));
        for (String roleDesc : roles) {
            Role role = subRoles.getByName(roleDesc);
            assertThat(roleDesc, role, not(nullValue()));
        }
        try (TabbedLineReader reader = new TabbedLineReader(initFile, 3)) {
            for (TabbedLineReader.Line line : reader) {
                String roleDesc = line.get(2);
                Role role = subRoles.getByName(roleDesc);
                assertThat(roleDesc, role.getId(), equalTo(line.get(0)));
            }
        }

    }

}
