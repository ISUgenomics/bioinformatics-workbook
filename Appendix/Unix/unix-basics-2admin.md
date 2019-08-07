---
title: "Administrative Commands"
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Unix Basics 2
This exercise will provide you details about some administrative commands with examples. Here, you can learn how to change permissions for files and folders to modify its accessibility and commands to obtain information about the system you are using.

## Changing permissions

All files in the UNIX system will have a set of permissions which define what can be done with that file and by whom. Here, what refers to read (view contents), write (modify) and execute (run as a script) and whom refers to user (owner), group (collection of users that the user belongs to) and others (everyone else).

| Permissions | Symbol    |
|:-------------|----------:|
| read         | `r`       |
| write        | `w`       |
| execute      | `x`       |
| all users    | `a`       |

| Relations | Symbol |
|:----------|-------:|
| owner     | `u`    |
| group     | `g`    |
| others    | `o`    |

To look at the permissions for any file, you can list the files with `l` option (`ls –l`).
Permissions	User	Group	Size	Date modified	Name

It looks something like this:
```
total 200
-rw-r--r--. 1 arnstrm domain users             174 Jan 15 23:36 1
drwxr-xr-x. 3 arnstrm domain users              38 Nov 15 10:31 Desktop
drwx------. 2 arnstrm domain users              10 Nov 15 10:38 Downloads
drwxr-xr-x. 3 arnstrm domain users              49 Jan  4 12:38 R
lrwxrwxrwx. 1 arnstrm domain users              17 Feb  3  2017 arnstrm -> /work/GIF/arnstrm
drwxr-xr-x. 3 arnstrm domain users             252 Mar 11  2017 bash_config
drwxr-xr-x. 2 arnstrm domain users              48 Aug 11 14:07 bin
drwxr-xr-x. 2 arnstrm domain users              10 Oct 25 12:15 ccp4_tmp
-rw-r--r--. 1 arnstrm domain users            4796 Jan 15 23:34 compnode
-rw-r--r--. 1 arnstrm domain users            3213 Jan 15 23:33 headnonde
-rwxr-xr-x. 1 arnstrm domain users          159656 Jan 15 23:47 ld-linux.so.2
-rw-r--r--. 1 arnstrm domain users             699 Oct 31 09:10 md5sum_severin
lrwxrwxrwx. 1 arnstrm domain users              12 Mar 20  2017 ncbi -> arnstrm/ncbi
-rw-r--r--. 1 arnstrm domain users            9968 Dec 31 14:35 ncbi_error_report.xml
-rw-------. 1 arnstrm domain users             522 Oct 30 14:45 nohup.out
-rw-r-----. 1 arnstrm domain users             287 Feb  7  2017 template.slurm
     -1-   -2-  -3-       -4-                  -5-       -6-       -7-
```
1. First letter of the first column specifies the type. It can be either `d` is directory, `l` is link or `-` is regular file. Remaining 9 letters of the first column, each 3 specifies permissions set for `user`, `group` and `others`, respectively. Here `r` is read, `w` is write, `x` is execute and `-` is blank or unset. The last `.` sign specifies attributes for this item (to see complete list go the official manual [here](https://www.gnu.org/software/coreutils/manual/html_node/What-information-is-listed.html))
2. Second column, specifies number of sub directories housed inside. It can also be number of links that points to it.
3. The owner of the file/directory: `user`
4. The fourth column `domain users` is the group, `user` belongs to.
5. Next, the number you see is the size of the listed entry. For a file, it will show the actualy size of the file in bytes, but for folder, it will not display them correctly as it wont consider the file sizes that are inside the directory.
6. The sixth column (`Jan 15 23:36`) is the month, day, and time on which the entry was last modified.
7. The last field, is the name of the listed entry.


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

### 8. IP address:

To get the IP address for the machine you can use the `ifconfig` command.

```
ifconfig
```
it will lists all properties as follows:

```
eth2      Link encap:Ethernet  HWaddr 00:0N:00:00:N0:NN
          inet addr:00.00.000.000  Bcast:00.00.000.000  Mask:000.000.000.0
          inet6 addr: 0000:000:000:000:00n:00nn:nn00:n000/00 Scope:Link
          UP BROADCAST RUNNING MULTICAST  MTU:9000  Metric:1
          RX packets:1693260659 errors:0 dropped:0 overruns:0 frame:0
          TX packets:600815878 errors:0 dropped:0 overruns:0 carrier:0
          collisions:0 txqueuelen:10000
          RX bytes:14812199265240 (13.4 TiB)  TX bytes:52487439229 (48.8 GiB)
...
<clipped rest of the output>
```
you can use the combination of commands to just display the IP address as follows:

```
ifconfig | grep -Eo 'inet (addr:)?([0-9]*\.){3}[0-9]*' | grep -Eo '([0-9]*\.){3}[0-9]*' | grep -v '127.0.0.1'
```


### 9. Other information:

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

---
[Table of contents](../programs.md)
