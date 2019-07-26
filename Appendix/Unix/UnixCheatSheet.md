---
title: "Unix CheatSheet"
layout: single
author: Arun Seetharam
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---


|Navigation| ||
|--|--|--|
| <span style="color:Green">_Command_</span>|<span style="color:Green">_Function_</span>|<span style="color:Green">_Syntax/example usage_</span> |
|`ls` 	|list contents	|`ls` <span style="color:Red">[OPTIONS] DIRECTORY</span>|
|`pwd`	|print working directory	|`pwd`|
|`cd`|change directory	|`cd` ~ or `cd` 		#home directory|
| | 			|`cd` .. #previous (parent directory)|

|File/Directory operations| | |
|--|--|--|
| <span style="color:Green">_Command_</span>|<span style="color:Green">_Function_</span>|<span style="color:Green">_Syntax/example usage_</span> |
|`mkdir`	|make directory	|`mkdir` <span style="color:Red">DIRECTORY</span>|
|`cp`|copy files/directories	|`cp `<span style="color:Red">SOURCE DESTINATION</span>|
|`man`|	manual page (help)|	`man `<span style="color:Red">COMMAND</span>|
|`mv`	|move files/directories	|`mv` <span style="color:Red">SOURCE DESTINATION</span>|
|`touch`|	create file	|`touch` <span style="color:Red">FILE</span>|
|`nano`	|edit file	|`nano` <span style="color:Red">FILE</span>|
|`less`	|view file (with more options)	|`less` <span style="color:Red">FILE</span>|
|`more`|view file (with less options)	|`more` <span style="color:Red">FILE</span>|
|`cat`	|catalog file contents	|`cat` <span style="color:Red">FILE</span>|
|`head`|show first few lines of a file	|`head` <span style="color:Red">FILE</span>|
|`tail`	|show last few lines of a file	|`tail` <span style="color:Red">FILE</span>|
|`rmdir`	|remove empty directory	|`rmdir` <span style="color:Red">DIRECTORY</span>|
|`rm`	|remove file(s)	|`rm` <span style="color:Red">FILE</span>|


|Compression/archiving| | |
|--|--|--|
| <span style="color:Green">_Command_</span>|<span style="color:Green">_Function_</span>|<span style="color:Green">_Syntax/example usage_</span> |
|`zip`	|zip compress	|`zip` <span style="color:Red">OUTFILE.zip INFILE.txt</span>|
|||`zip` `-r` <span style="color:Red">OUTDIR.zip DIRECTORY</span>|
|`unzip`	|decompress zipped file	|`unzip` <span style="color:Red">ANYTHING.zip</span>|
|`tar`	|archive and compress files/directories	|`tar` `-czvf` <span style="color:Red">OUTFILE.tar.gz DIRECTORY</span> #compress|
|||`tar` `-xzvf` <span style="color:Red">OUTFILE.tar.gz</span> # extract|
|`gzip`	|gzip files	|`gzip` <span style="color:Red">SOMEFILE</span>|
|`gunzip`	|decompress gzipped files	|`gunzip` <span style="color:Red">SOMEFILE.gz</span>|

|File permissions| | |
|--|--|--|
| <span style="color:Green">_Command_</span>|<span style="color:Green">_Function_</span>|<span style="color:Green">_Syntax/example usage_</span> |
|`chmod`	|change permissions for files/directories	|`chmod` <span style="color:Red">[OPTIONS] RELATIONS[+/-]PERMISSIONS FILE</span>|


|File manipulations| | |
|--|--|--|
| <span style="color:Green">_Command_</span>|<span style="color:Green">_Function_</span>|<span style="color:Green">_Syntax/example usage_</span> |
|`grep`	|search a pattern	|grep <span style="color:Red">[OPTIONS] "PATTERN" FILENAME</span>|
|`sed`	|stream edit a file	|`sed 's/search/replace/g'` <span style="color:Red">FILENAME</span>|
|`awk`	|multi-purpose command	|`awk 'PATTERN {ACTION}'` <span style="color:Red">FILENAME</span>|
|`tr`	|translate or transliterate a file	|`tr` <span style="color:Red">[OPTIONS] "STRING1" "STRING2"</span>  `<` <span style="color:Red">INFILE</span> |
|`wc`	|word count	|`wc` <span style="color:Red">FILENAME</span>|
|`sort`	|sort files	|`sort` <span style="color:Red">FILE1<\span> `>` <span style="color:Red">SORTED_FILE1</span>|
|`uniq`	|display unique lines 	|`uniq` <span style="color:Red">[OPTIONS] INFILE</span> `>` <span style="color:Red">OUTFILE</span>|
|`diff`	|display difference	|`diff` <span style="color:Red">[OPTIONS] FILE1 FILE2</span>|
|`comm`	|display common lines among files	|`comm` <span style="color:Red">[OPTIONS] FILE1 FILE2</span>|
|`cut`	|break files vertically based on fields	| `cut` `–d` <span style="color:Red">"DELIMITER"</span> `–f` <span style="color:Red">NUMBER FILE</span>|
|`split`	|break files horizontally 	|split <span style="color:Red">[OPTIONS] FILENAME</span>|
|`paste`	|combine files side by side	|`paste` <span style="color:Red">FILE1 FILE2</span> `>` <span style="color:Red">FILE3</span>|
|`join`	|join files based on common field	|`join -t` <span style="color:Red">'DELIMITER'</span> `-1` <span style="color:Red">N</span> `-2` <span style="color:Red">N FILE1 FILE2</span>|
