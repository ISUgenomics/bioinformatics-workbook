# Unix Basics 2
This exercise will provide you details about some administrative commands with examples. Here, you can learn how to change permissions for files and folders to modify its accessibility and commands to obtain information about the system you are using.

## Changing permissions

All files in the UNIX system will have a set of permissions which define what can be done with that file and by whom. Here, what refers to read (view contents), write (modify) and execute (run as a script) and whom refers to user (owner), group (collection of users that the user belongs to) and others (everyone else). 

| Permissions | Letter    |
|:-------------|----------:|
| read         | `r`       |
| write        | `w`       |
| execute      | `x`       |
| all users    | `a`       |

| Relations | Letter |
|:----------|-------:|
| owner     | `u`    |
| group     | `g`    |
| others    | `o`    |

To look at the permissions for any file, you can list the files with `l` option (`ls –l`). 
Permissions	User	Group	Size	Date modified	Name
 
assets/homedir.png 

(d=directory, l=link, r=read, w=write, x=execute, -=blank, u=user, g=group, o=others)

To set/modify a file's permissions you need to use the `chmod` command (`ch`ange `mod`e). Only the owner of a file can alter a file's permissions. The syntax: 
```
chmod [OPTIONS] RELATIONS[+ or -]PERMISSIONS FILE
```
### 1. Adding permissions

```
chmod RELATIONS+PERMISSIONS FILENAME
```
is the add permissions syntax, an example
```
chmod g+rwx FILENAME
```
which grants read, write and execute permissions for group

```
chmod g+r FILENAME
```
grants read permission for group

```
chmod a+rwx FILENAME
```
makes the file public (don’t do this to any file/directory unless you want to share)


### 2. Removing permissions:

```
chmod RELATIONS-PERMISSIONS FILENAME
```
is the syntax for removing permissions, examples:

```
chmod g-wx FILENAME
```
removes write and execute permissions for group
```
chmod g-rwx FILENAME
```
removes all permissions for group
```
chmod a-rwx FILENAME
```
removes all permissions for others
```
chmod a-x FILENAME
```
removes execution permissions for others

OPTIONS include `-R` recursively (the permissions are applied to all the files, directories present inside the directory)


***Task 1.12: Check the permissions for the files located in the tutorials directory.***
```
ls -l
```
*What permissions does the group have on these files? Which group does your account belong to?*

## Check system properties

In this section, you will learn how to check system resources (space, memory, disk usage, storage properties), system properties (operating system, Linux version, kernel version) and commands to access other information (CPU type, memory type, variables available etc ) about the environment

### 1. Directory size

To get the size of the directory, you can use the `du` command (`d`isk `u`sage)
```
du -sh DIRECTORY
```
The options are to summarize (`s`) and human readable format (`h`). While the summarize will avoid printing size for every file in the directory, the human readable format will give the folder size in kilo/mega/giga bytes  (instead of bytes).

### 2. File size

If you are interested in knowing the size of a particular file, you can use the `ls` command with `l` and `h` options 

```
ls -lh FILENAME
```
Check unix-basics-1 for more details about list command. The option `l` list in a long format and `h` in human readable file sizes.

### 3. Available storage and mounts

To display free disk space and mounted devices, you can use the `df` command. If no file name is given, the space available on all currently mounted file systems is shown

```
df -h
```
again, using the option `-h` will give you results in human readable format.


### 4. Available memory

If you want to see how much memory is available on your machine, you can use the `free` command.
```
free
```
Output would be:
```
             total       used       free     shared    buffers     cached
Mem:      32874744   32607664     267080          0      77600   31013192
-/+ buffers/cache:    1516872   31357872
Swap:     61438900     873856   60565044
```
As you can see the numbers are in bytes and very difficult to understand. You can modify this default behavior using some options. some options to modify this are:


| Options | What it does                         |
|--------:|:-------------------------------------|
| `-g`    | display numbers in gigabytes         |
| `-m`    | display numbers in megabytes         |
| `-k`    | display numbers in kilobytes         |
| `-b`    | display numbers in bytes, (default). |

Example:
```
free -g
```
output would be:
```
             total       used       free     shared    buffers     cached
Mem:            31         31          0          0          0         29
-/+ buffers/cache:          1         29
Swap:           58          0         57
```
This is much easier to understand.

### 5. System properties


Just to get the Operating system name:
```
cat /etc/system-release
```
Output would be:
```
Red Hat Enterprise Linux Server release 6.4 (Santiago)
```

If no such file, they try:
```
cat /etc/*release*
```
you might get:
```
Red Hat Enterprise Linux Server release 6.4 (Santiago)
Red Hat Enterprise Linux Server release 6.4 (Santiago)
cpe:/o:redhat:enterprise_linux:6server:ga:server
```

The other command is the `uname`
```
uname -a 
```
output would be:
```
Linux hpc5 2.6.32-358.11.1.el6.x86_64 #1 SMP Wed May 15 10:48:38 EDT 2013 x86_64 x86_64 x86_64 GNU/Linux
```
which is Kernel, node, kernel version, kernel release date, machine type, processor type, platform and OS type, respectively.

You can also ask for a specific thing by using these options:

```
uname -s # kernel name
uname -n # node
uname -v # version
uname -r # release version date
uname -i # platform
uname -m # machine type
uname -p # processor type
uname -o # OS type
```

### 7. Processor and Memory information:

These information will be in the file. Just by cataloging the file, you can find read these information:
```
cat /proc/meminfo
```
for the memory information, and

```
cat /proc/cpuinfo
```
for the CPU information

### 8. Other information:

For getting more information about the environment, you can type `env`, which lists all the variables currently set. If you want to know specifically about a variable, you can do:

```
echo $VARIABLE
```
Some variables that are useful are:

| Variable   | Information                        |
|:-----------|:-----------------------------------|
| `HOSTNAME` | hostname for the system            |
| `TERM`     | terminal                           |
| `SHELL`    | Shell type (bash, csh, ksh etc)    |
| `USER`     | Username                           |
| `PATH`     | paths where executables are stored |
| `PWD`      | present working directory          |
| `EDITOR`   | default text editor                |
| `HOME`     | path for home                      |
| `DISPLAY`  | where to route the display         |
| `HISTFILE` | file where the history is saved    |

