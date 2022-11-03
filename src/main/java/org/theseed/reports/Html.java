/**
 *
 */
package org.theseed.reports;

import static j2html.TagCreator.*;

import java.util.Collection;
import java.util.List;

import j2html.tags.ContainerTag;
import j2html.tags.DomContent;

/**
 * These are static methods to help when building web pages.
 *
 * @author Bruce Parrello
 *
 */
public class Html {

    // STYLE CONSTANTS

    /** background color for values indicating bad genomes */
    public static final String BAD_STYLE = "background-color: gold;";

    public static final String CSS_HREF = "https://core.theseed.org/SEEDtk/css/p3.css";

    private static final String NUM_FORMAT = "%-8.2f";

    private static final String INT_FORMAT = "%d";

    public static final String TABLE_CLASS = "p3basic";

    public static final String BODY_CLASS = "claro";

    public static final String EXTRA_STYLES = "	table." + TABLE_CLASS + " th {\n" +
            "		background-color: #c8c8c8;\n" +
            "	}\n" +
            "	table." + TABLE_CLASS + " th.num, table." + TABLE_CLASS + " td.num {\n" +
            "		text-align: right;\n" +
            "	}\n" +
            "	table." + TABLE_CLASS + " th.flag, table." + TABLE_CLASS + " td.flag {\n" +
            "		text-align: center;\n" +
            "	}\n" +
            "   h1, h2 {\n" +
            "		font-weight: bolder;\n" +
            "	}\n" +
            "   h1, h2, p, ul {\n" +
            "		margin: 12px 12px 0px 12px;\n" +
            "	}\n" +
            "   table.p3basic ul {\n" +
            "		margin: 0px;\n" +
            "		list-style: disc outside none;\n" +
            "		padding-left: 20px;\n" +
            "	}\n" +
            "   table.p3basic li {\n" +
            "		margin: 3px 0px;\n" +
            "	}\n" +
            "   div.wrapper {\n" +
            "       margin: 12px;\n" +
            "   }\n" +
            "   div.shrinker {\n" +
            "       margin: 12px;\n" +
            "       display: inline-block;\n" +
            "       min-width: 0;\n" +
            "       width: auto;\n" +
            "   }\n" +
            "   li {\n" +
            "		margin: 6px 12px 0px 12px;\n" +
            "	}\n" +
            "	table." + TABLE_CLASS + " {\n" +
            "	    display:table;\n" +
            "	}\n";


    // LINK CONSTANTS

    private static final String PAGE_LINK_FMT = "%s.html";

    /**
     * Display the specified value in an alternate color if the flag is FALSE.
     *
     * @param flag	TRUE if the value is good, else FALSE
     * @param str	text to display
     */
    public static ContainerTag colorCell(boolean flag, String str) {
        DomContent cell;
        if (str.isEmpty())
            cell = rawHtml("&nbsp;");
        else
            cell = text(str);
        ContainerTag retVal = colorCell(flag, cell);
        return retVal;
    }

    /**
     * Display an empty table cell.
     */
    public static ContainerTag emptyCell() {
        return td(rawHtml("&nbsp;"));
    }

    /**
     * Display the specified HTML in an alternate color if the flag is FALSE.
     *
     * @param flag	TRUE if the value is good, else FALSE
     * @param cell	HTML to display
     */
    public static ContainerTag colorCell(boolean flag, DomContent cell) {
        ContainerTag retVal = td(cell);
        if (! flag) {
            retVal = retVal.withStyle(BAD_STYLE);
        }
        return retVal;
    }

    /**
     * @return a link to the specified genome's report page
     *
     * @param genome_id		ID of the target genome
     * @param text			text for the link
     */
    public static DomContent gPageLink(String genome_id, String text) {
        return a(text).withHref(String.format(PAGE_LINK_FMT, genome_id));
    }

    /**
     * Display the specified value in an alternate color if the flag is FALSE.
     *
     * @param flag	TRUE if the value is good, else FALSE
     * @param val	floating-point value to display
     */
    public static ContainerTag colorCell(boolean flag, double val) {
        return colorCell(flag, num(val)).withClass("num");
    }

    /**
     * Display the specified value in an alternate color if the flag is FALSE.
     *
     * @param flag	TRUE if the value is good, else FALSE
     * @param val	integer value to display
     */
    public static ContainerTag colorCell(boolean flag, int val) {
        return colorCell(flag, num(val)).withClass("num");
    }

    /**
     * Display the specified floating-point value in a table cell.
     *
     * @param val	floating-point value to display
     */
    public static ContainerTag numCell(double val) {
        return td(num(val)).withClass("num");
    }

    /**
     * Display the specified integer value in a table cell.
     *
     * @param val	floating-point value to display
     */
    public static ContainerTag numCell(int val) {
        return td(num(val)).withClass("num");
    }

    /**
     * Display a flag value in a cell.  A false value is in the alternate color.
     *
     * @param flag		TRUE if the relevant attribute is acceptable, else FALSE
     * @param trueText	text to display if the flag is TRUE
     * @param falseText	test to display if the flag is FALSE
     */
    public static ContainerTag flagCell(boolean flag, String trueText, String falseText) {
        String text = (flag ? trueText : falseText);
        ContainerTag retVal = colorCell(flag, text);
        return retVal.withClass("flag");
    }

    /**
     * Create a web page with the specified title and body components.
     *
     * @param pTitle	title for the page
     * @param bodyItems	one or more items for the body of the page
     *
     * @return the HTML string for the web page
     */
    public static String page(String pTitle, DomContent... bodyItems) {
        String retVal = html(
                head(
                        title(pTitle),
                        link().withRel("stylesheet").withHref(CSS_HREF).withType("text/css"),
                        style(EXTRA_STYLES).withType("text/css")
                    ),
                    body(bodyItems).withClass(BODY_CLASS)
                ).render();
        return retVal;
    }

    /**
     * Add a row to the statistical details table row collection.
     *
     * @param detailRows	statistical details table row collection
     * @param label		label for the row
     * @param cell			table cell with the data
     */
    public static void detailRow(List<DomContent> detailRows, String label, ContainerTag cell) {
        detailRows.add(tr(th(label), cell));
    }

    /**
     * @return a formatted integer
     *
     * @param val	integer to format
     */
    public static String num(int val) {
        return String.format(INT_FORMAT, val);
    }

    /**
     * @return a formatted floating-point number
     *
     * @param val	floating-point number to format
     */
    public static String num(double val) {
        return String.format(NUM_FORMAT, val);
    }

    /**
     * @return a table created from the specified detail rows
     *
     * @param collection		rows to put in the table
     */
    public static DomContent formatTable(String header, Collection<DomContent> collection) {
        return join(
                h2(header),
                div(table().with(collection.stream()).withClass(TABLE_CLASS + " striped")).withClass("wrapper")
            );
    }

}
