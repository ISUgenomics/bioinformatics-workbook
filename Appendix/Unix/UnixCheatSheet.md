---
title: "Unix CheatSheet"
layout: single
author: Arun Seetharam
author1: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

## Navigation

|Navigation| ||
|--|--|--|
| <span style="color:Green">_Command_</span>|<span style="color:Green">_Function_</span>|<span style="color:Green">_Syntax/example usage_</span> |
|`ls` 	|list contents	|`ls` <span style="color:Red">[OPTIONS] DIRECTORY</span>|
|`pwd`	|print working directory	|`pwd`|
|`cd`|change directory	|`cd` ~ or `cd` 		#home directory|
| | 			|`cd` .. #previous (parent directory)|

## File/Directory operations

|File/Directory operations | | |
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


## File manipulations

|File manipulations| | |
|--|--|--|
| <span style="color:Green">_Command_</span>|<span style="color:Green">_Function_</span>|<span style="color:Green">_Syntax/example usage_</span> |
|`grep`	|search a pattern	|grep <span style="color:Red">[OPTIONS] "PATTERN" FILENAME</span>|
|`sed`	|stream edit a file	|`sed 's/search/replace/g'` <span style="color:Red">FILENAME</span>|
|`awk`	|multi-purpose command	|`awk 'PATTERN {ACTION}'` <span style="color:Red">FILENAME</span>|
|`tr`	|translate or transliterate a file	|`tr` <span style="color:Red">[OPTIONS] "STRING1" "STRING2"</span>  `<` <span style="color:Red">INFILE</span> |
|`wc`	|word count	|`wc` <span style="color:Red">FILENAME</span>|
|`sort`	|sort files	|`sort` <span style="color:Red">FILE1</span> `>` <span style="color:Red">SORTED_FILE1</span>|
|`uniq`	|display unique lines 	|`uniq` <span style="color:Red">[OPTIONS] INFILE</span> `>` <span style="color:Red">OUTFILE</span>|
|`diff`	|display difference	|`diff` <span style="color:Red">[OPTIONS] FILE1 FILE2</span>|
|`comm`	|display common lines among files	|`comm` <span style="color:Red">[OPTIONS] FILE1 FILE2</span>|
|`cut`	|break files vertically based on fields	| `cut` `–d` <span style="color:Red">"DELIMITER"</span> `–f` <span style="color:Red">NUMBER FILE</span>|
|`split`	|break files horizontally 	|split <span style="color:Red">[OPTIONS] FILENAME</span>|
|`paste`	|combine files side by side	|`paste` <span style="color:Red">FILE1 FILE2</span> `>` <span style="color:Red">FILE3</span>|
|`join`	|join files based on common field	|`join -t` <span style="color:Red">'DELIMITER'</span> `-1` <span style="color:Red">N</span> `-2` <span style="color:Red">N FILE1 FILE2</span>|

## Compression/archiving

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

## File permissions

|File permissions| | |
|--|--|--|
| <span style="color:Green">_Command_</span>|<span style="color:Green">_Function_</span>|<span style="color:Green">_Syntax/example usage_</span> |
|`chmod`	|change permissions for files/directories	|`chmod` <span style="color:Red">[OPTIONS] RELATIONS[+/-]PERMISSIONS FILE</span>|


## ADDITIONAL COMMANDS

|ADDITIONAL COMMANDS| |
|--|--|
| <span style="color:Green">_Command_</span>|<span style="color:Green">_Function_</span>|
|`du –sh` |DIR	show directory size|
|`whoami`	|display username|
|`date`	|system date/time|
|`cal`	|calendar|
|`find . –name FILE`	|find a file/directory|
|`which CMD`	|display default cmd path|
|`whereis CMD`	|show possible locations of cmd|
|`locate FILE`	|find instances of a file|
|`clear`	|clear screen|
|`sleep 5`	|pause 5 (any) seconds|
|`top`	|current running processes |
|`ps`	|current running processes|
|`wget` |URL	download specified URL|

## SHORTCUTS

|SHORTCUTS| |
|--|--|
| <span style="color:Green">_Command_</span>|<span style="color:Green">_Function_</span>|
|`TAB`	|autocomplete names|
|`UP/DOWN`	|browse previous commands|
|`ctrl+c`	|interrupt/kill anything|
|`ctrl+l`	|clear screen|
|`ctrl+d`	|quit, exit|
|`ctrl+z`	|suspend (use fg to restore)|
|`!!`	|repeat last command|
|`alt+.`	|last argument of previous cmd|
|`ctrl+insert`	|copy selection|
|`shift+insert`	|paste copied text|
|`ctrl+a`	|go to start of the line|
|`ctrl+e`	|go to end of the line|
|`ctrl+r`	|reverse search history|
|`cd ~`	|go to home|

## NANO SHORTCUTS

|NANO SHORTCUTS| |
|--|--|
| <span style="color:Green">_Command_</span>|<span style="color:Green">_Function_</span>|
|`ctrl+r`	|read/insert file|
|`ctrl+o`	|save file|
|`ctrl+x`	|close file|
|`alt+a`	|start selecting text|
|`ctrl+k`	|cut selection|
|`ctrl+u`	|uncut (paste) selection|
|`alt+/`	|go to end of the file|
|`ctrl+a`	|go to start of the line|
|`ctrl+e`	|go to end of the line|
|`ctrl+c`	|show line number|
|`ctrl+_`	|go to line number |
|`ctrl+w`	|find matching word|
|`alt+w`	|find next match|
|`ctrl+\`	|find and replace|

## PIPES, REDIRECTS

|PIPES, REDIRECTS| |
|--|--|
| <span style="color:Green">_Command_</span>|<span style="color:Green">_Function_</span>|
|cmd `<` <span style="color:Red">file</span>	|use file as input|
|cmd `>` <span style="color:Red">file</span>	|write output to file|
|cmd `>>` <span style="color:Red">file</span>	|append output to file|
|cmd `2>` <span style="color:Red">stderr</span>	|error output to file|
|cmd `1>&2` <span style="color:Red">file</span>	|send output and error to file|
|cmd1 `\|` <span style="color:Red">cmd2</span> 	|send output of cmd1 to cmd2|

## PRE-DECLARED VARIABLES

|PRE-DECLARED VARIABLES| |
|--|--|
|<span style="color:Green">Variables\*</span>	|<span style="color:Green">Description</span>|
|`$USER`	|username|
|`$HOME`|	home path|
|`$PWD`|	working directory path|
|`$PATH`|	path for executables|
|`$HOSTNAME`|	machine name|
|`$SHELL`|	current shell|
|`$SSH_CLIENT`|	local client's IP address|
|`$TERM`	|type of terminal|

* env command lists all the assigned variables

## HPC-CLUSTER Commands for SLURM

|HPC-CLUSTER Commands | |TORQUE|
|--|--|--|
|<span style="color:Green">Command</span>	|<span style="color:Green">Function</span>	|<span style="color:Green">Syntax/example usage</span> |
|sinfo -a	|list all queues|sinfo -a|
|squeue	|list all jobs|squeue|
|squeue -u userid	|list jobs for userid|squeue -u userid|
|squeue -t R	|list running jobs|squeue -t R|
|smap	|show jobs, partitions and nodes in a graphical network topology|smap|
|sbatch	|submit a slurm job	|sbatch [script]|
|scancel|	delete slurm batch job|	scancel [job_id]	|
|scontrol hold	|hold slurm batch jobs	|scontrol hold [job_id]|
|scontrol release	|release hold on slurm batch jobs	|scontrol release [job_id]	|

* [Tutorial on SLURM commands](../HPC/SLURM/slurm-cheatsheat.md)

## HPC-CLUSTER Commands for TORQUE

|HPC-CLUSTER Commands | |TORQUE|
|--|--|--|
|<span style="color:Green">Command</span>	|<span style="color:Green">Function</span>	|<span style="color:Green">Syntax/example usage</span> |
|`squeue`	|show state of jobs	|`squeue` `–a`		# current jobs on cluster |
|`squeue` `–u` <span style="color:Red">username</span>	||# current jobs by the user|
|`squeue` `–j` <span style="color:Red">jobid</span>	||# information about the job (id#)|
|`squeue` `-l –u` <span style="color:Red">username</span>	||# current jobs by the user|
|`scancel`	|delete job from the queue	|`scancel` <span style="color:Red">jobid</span>|
|`sbatch`	|submit job to the queue	|`sbatch` <span style="color:Red">submissionfile.sub</span>
|`scontrol`	|control jobs	|`scontrol hold` <span style="color:Red">jobid jobid</span> # hold the job
|||`scontrol release` <span style="color:Red">jobid jobid</span> # release the job|
|||`scontrol show` <span style="color:Red">jobid jobid</span> # info on the job|
|`srun`	|run a job command	|`srun` `-N 1 -n 16 -t 4:00:00 --pty bash` # start a interactive job session|
|`sinfo`|	show state of nodes and partitions	|`sinfo`|
|||`sinfo` `-p` <span style="color:Red">tiny</span>      # show info on tiny partition|
|`smap`	|show state of jobs, nodes and partitions (colored) 	|`smap`|
|`module`	|use preinstalled programs	|`module load` <span style="color:Red">PROGRAM</span>		# loads program for use|
|||`module list`			# lists all loaded modules|
|||`module avail`			# lists available modules|
|||`module unload` <span style="color:Red">PROGRAM</span>	# unloads module|

PS: An A-Z Index of the Bash command line for Linux can be found at found [Here](http://ss64.com/bash/index.html)
