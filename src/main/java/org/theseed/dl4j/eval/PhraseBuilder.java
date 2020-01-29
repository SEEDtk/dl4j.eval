/**
 *
 */
package org.theseed.dl4j.eval;

import java.util.ArrayList;

/**
 * This is a utility class that allows us to build a phrase from components, inserting proper punctuation.
 * The constructor specifies a prefix and a suffix.  If no phrases are added, only the suffix is returned.
 * If any phrases are added, the prefix is put in first.  The number of phrases determines how they are
 * connected.
 *
 * @author Bruce Parrello
 */

public class PhraseBuilder {

    // FIELDS

    /** list of phrases */
    private ArrayList<String> phrases;
    /** prefix to the first phrase */
    private String prefix;
    /** suffix for the end */
    private String suffix;

    /**
     * Construct an empty modifier phrase.
     *
     * @param prefix	prefix to be used if at least one phrase is added
     * @param suffix	suffix to put at the end
     */
    public PhraseBuilder(String prefix, String suffix) {
        this.phrases = new ArrayList<String>();
        this.prefix = prefix;
        this.suffix = suffix;
    }

    /**
     * Accumulate a new phrase.
     *
     * @param phrase	phrase to add
     */
    public void addPhrase(String phrase) {
        this.phrases.add(phrase);
    }

    /**
     * Assemble and return the modifier phrase.
     */
    @Override
    public String toString() {
        StringBuilder retVal = new StringBuilder();
        int phraseCount = this.phrases.size();
        if (phraseCount > 0) {
            // We are non-empty. Start with the prefix.
            retVal.append(this.prefix);
            retVal.append(' ');
            switch (phraseCount) {
            case 2:
                // Two phrases are joined with " and ".
                retVal.append(this.phrases.get(0));
                retVal.append(" and ");
                retVal.append(this.phrases.get(1));
                break;
            case 1:
                // A single phrase is unadorned.
                retVal.append(this.phrases.get(0));
                break;
            default:
                // Three or more requires an Oxford comma.
                retVal.append(this.phrases.get(0));
                int lastPhrase = phraseCount - 1;
                for (int i = 1; i < lastPhrase; i++) {
                    retVal.append(", ");
                    retVal.append(this.phrases.get(i));
                }
                retVal.append(", and ");
                retVal.append(this.phrases.get(lastPhrase));
            }
        }
        retVal.append(this.suffix);
        return retVal.toString();
    }



}
