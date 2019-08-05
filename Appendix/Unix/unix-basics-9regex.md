# Regular Expressions


Imagine, rather than typing all the file names, you could represent them with a single word.  This word or a pattern used to represent files/directories are called regular expression (regex for short).
* **Simple examples:**
*  to represent any word: eg. \*.txt means all files with txt extension
  * ? to represent a single letter: eg. ?????.txt matches all files with exactly 5 letters, with txt extension.
  * ^ beginning and $ for the end of the word:  eg  ^text* forces the match for the beginning letters only.
