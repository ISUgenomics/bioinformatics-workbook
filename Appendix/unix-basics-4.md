---
title: "Useful Programs and Unix Basics"
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

## All about SED command

The  `s`treamline `ed`itor  or `sed` command is a stream editor that reads one or more text files, makes changes or edits according to editing script, and writes the results to standard output. First, we will discuss `sed` command with respect to search and replace function. Other uses for the `sed` can also be found in [this](http://www.grymoire.com/Unix/Sed.html#uh-47) official guide, but we will discuss them briefly in this chapter as well.

### 1. Find and replace
Most common use of `sed` is to substitute text, matching a pattern. The syntax for doing this in `sed` is as follows:

```
sed 'OPERATION/REGEXP/REPLACEMENT/FLAGS' FILENAME
```

  - Here, `/` is the delimiter (you can also use `_` (underscore), `|` (pipe) or `:` (colon) as delimiter as well)
  - `OPERATION` specifies the action to be performed (sometimes if a condition is satisfied). The most common and widely used operation is `s` which does the substitution operation (other useful operators include `y` for transformation, `i` for insertion, `d` for deletion etc.).
  - `REGEXP` and `REPLACEMENT` specify search term and the substitution term respectively for the operation that is being performed.
  - `FLAGS` are additional parameters that control the operation. Some common `FLAGS` include:
      * `g`	replace all the instances of `REGEXP` with `REPLACEMENT` (globally)
      * `N` where N is any number, to replace Nth instance of the `REGEXP` with `REPLACEMENT`
      * `p` if substitution was made, then prints the new pattern space
      * `i` ignores case for matching `REGEXP`
      * `w` file If substitution was made, write out the result to the given file
      * `d` when specified without `REPLACEMENT`, deletes the found `REGEXP`

Some search and replace examples:

find and replace all chr to chromosome in the file
```
sed 's/chr/chromosome/g' FILENAME > NEWFILE
```
find and replace, but only the one instance per line (first occurrence of chr will be changed to chromosome)
```
sed 's/chr/chromosome/1' FILENAME > NEWFILE
```
find and replace, but do it directly on the original file
```
sed -i 's/chr/chromosome/g' FILENAME
```
find and replace directly, but save a old version too
```
sed -i.old 's/chr/chromosome/g' FILENAME
```
find and replace, only if you also find MTF1 in the line
```
sed '/MTF1/s/chr/chromosome/g' FILENAME > NEWFILE
```

### 2. Print specific lines of the file

To print a specific line, you can use the address function, note that by deafault, `sed` will stream the entire file, so when you are interested in specific lines only, you will have to suppress this feature using the option `-n`.

print 10th line
```
sed -n '10p' FILENAME
```
You can provide any number of additional lines to print using `-e` option (you can add any number of lines like this)

print 10th and 15th line
```
sed -n -e '10p' -e '15p' FILENAME
```
It also accepts range, using `,`

print lines 10 to 50
```
sed -n '10,50p' FILENAME
```
or you can create specific pattern, like multiple of a number using `~`

Every tenth line starting from 10, 20, 30.. to end of the file
```
sed -n '10~10p' FILENAME
```
print odd-numbered lines
```
sed -n '1~2p' FILENAME
```

Most powerful feature is that you can combine these ranges or multiples in any fashion. Example: `fastq` files have header on first line and sequence in second, next two lines will have the quality and a blank extra line (four lines make one read). Sometimes you will only need the sequence and header

to print 1,2,5,6,9,10 so on you can use
```
sed -n '1~4p;2~4p' FASTQ_FILE
```
pipe this to make a fasta file
```
sed -n '1~4p;2~4p' FASTQ_FILE | sed 's/^@/>/g' > FASTA_FILE
```
More combinations:

print 1 to 10, and then multiples of 10
```
sed -n '1,10~10p'  FILENAME
```

### 3. Delete specific lines of the file

All the above address types (specific line, range, multiples), also works with other types of operation, such as deletion and insertion. For deletion, you need to swap `p` with `d`

delete first line
```
sed "1d" FILENAME
```
delete lines 1 thru 3
```
sed "1,3d" FILENAME
```
delete blank lines
```
sed 's/^$//g' FILENAME
```

### 4. Insert specific lines to a file
Here, you use `i` for inserting text anywhere in the file


put "line to insert" in the second line
```
sed '2 i line to insert' FILENAME
```

### 5. Extract portions of the lines
still to come!

---
[Table of contents](programs.md)
