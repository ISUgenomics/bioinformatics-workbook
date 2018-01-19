## Advanced Commands

### grep

`grep` (`g`lobally search a `r`egular `e`xpression and `p`rint) is one of the most useful commands in UNIX and it is commonly used to filter a file/input, line by line, against a pattern eg., to print each line of a file which contains a match for pattern.
```
grep PATTERN FILENAME
```
Like any other command there are various options available for this command. Most useful options include:

| Argument  | Function                                                          |
|----------:|:------------------------------------------------------------------|
| `-v`      | inverts the match or finds lines NOT containing the pattern.      |
| `--color` | colors the matched text for easy visualization                    |
| `-F`      | interprets the pattern as a literal string.                       |
| `-H`,`-h` | print, don't print the matched filename                           |
| `-i`      | ignore case for the pattern matching.                             |
| `-l`      | lists the file names containing the pattern (instead of match).   |
| `-n`      | prints the line number containing the pattern (instead of match). |
| `-c`      | counts the number of matches for a pattern                        |
| `-o`      | only print the mattching pattern 			                        |
| `-w`      | forces the pattern to match an entire word.                       |
| `-x`      | forces patterns to match the whole line.                          |

With options, syntax is

```
grep [OPTIONS] PATTERN FILENAME
```

Some typical scenarios to use grep:
  -	Counting number of sequences in a multi-fasta sequence file
  -	Get the header lines of fasta sequence file
  -	Find a matching motif in a sequence file
  -	Find restriction sites in sequence(s)
  -	Get all the Gene IDs from a multi-fasta sequence files and many more.

Now let's use `grep` command to do some simple jobs with the sequences:

1. *Counting sequences*: By `FASTA` format definition, we know that number of sequences in a file should be equal to the number of description lines. So by counting `>` in file, you can count the number of sequences. This can be done using counting option of the `grep` with its count option `-c`.
```
grep -c ">" FILENAME
```
However, note that if the deflines somehow have `>` more than once, it will mess up the count! to be safe, you can use:
```
grep -c "^>" FILENAME
```
***Task 2.3: Count the number of sequences AT_cDNA.fa and RefSeq.faa***
```
grep -c ">" AT_cDNA.fa
grep -c ">" RefSeq.faa
```

2. *Looking for information*: If you are looking for information about the sequences, you can list all the headers (description lines) for the sequences using grep. Simply search for `>` and grep will list all the description lines.
```
grep ">" FILENAME
grep ">" AT_cDNA.fa
```
Alternatively, you can send it to a file if you want to use it later or you can just pipe it to less or more command to scroll through it line by line or page by page.
```
grep ">" FILENAME > HEADERFILE.txt
grep ">" FILENAME | less
grep ">" AT_cDNA.fa | less		
```
Use  &#8593; or &#8595; arrow keys to move up and down, press `q` to exit

3. *Subtracting two list of gene ids:* If there is a small list of genes that you want to remove from a larger list, you can use the grep function with these options:
```
grep -Fvw -f sub_list.txt full_list.txt
```
here `-F` and `-w` will make sure that the full word is used as literal string, `-v` will NOT print the matching patterns and `-f filename.txt` is to say that the input patterns are in the file.

4. *Count a word:* Unlink previous example, if the word your are searching occurs more than once in a line, it will only be counted once. To avoid this, you need to use a special option
```
grep -o "PATTERN" FILENAME
```

Now, let us have some fun with `grep`! See what kind of sequences are in `AT_cDNA.fa` file. Do they all seem to belong to same organism? Which organism?

Using `grep` you can also locate all the lines that contain a specific term you are looking for. This is very useful, especially to look for a specific gene among a large number of annotated sequences.
```
grep "word or phrase to search" FILENAME
```
***Task 2.4: Try searching for your favorite gene, to see if it is present in AT_cDNA.fa (this file contains all annotated sequences for Arabidopsis thaliana). Unlike Google or any search engines, only exact search terms will be identified, but you can ask grep to ignore cases while searching using -i option. Try these:***
```
grep -i "transcription factor" AT_cDNA.fa
grep -i "TFIIIA" AT_cDNA.fa
```
You can also use this feature to see if your sequence of interest has a specific feature (restriction site, motif etc.,) or not. This can be performed better using `--color` option of the `grep`.

Go to the sequences directory, search for `EcoR1` (`GAATTC`) site in the `NT21.fa` file, and use the color option. Also, try looking for a `C2H2` zinc finger motif in `RefSeq.faa` file (for simplicity let's assume zinc finger motif to be `CXXXCXXXXXXXXXXHXXXH`. Either you can use dots to represent any amino acids or use complex regular expressions to come up with a more representative pattern. Try these:
````
grep --color "GAATTC" ./Sequences/NT21.fa
grep --color "C..C............H...H" RefSeq.faa
```

You can also use `grep` command to exclude the results containing your search term. Say if you want to look at genes that are not located in chromosome 1, you can exclude it form your search by specifying `-v` option.
```
grep -i "transcription factor" AT_cDNA.fa| grep -v "chr1"
grep -i "transcription factor" AT_cDNA.fa| grep "chr1"
```
Notice the difference in output from the above two commands.
Try to understand the following command lines (and record your results, where applicable):
```
grep -c -w "ATP" RefSeq.faa
grep -c CGT[CA]GTG AT_cDNA.fa
grep -l "ATG" ./sequences/*.fa
```
You can also try some regular expressions related to nucleotide/protein sequences provided earlier to see how it works.
