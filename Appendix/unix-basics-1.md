# Unix Basics 1
This exercise is designed to provide the basic skills required for working in the UNIX environment, using plenty of relevant examples, specifically for biologists.  If you are using your personal computer, make sure that you have downloaded the files required for the workshop. This exercise will provide you information regarding navigation, files and directory creation/modification and some administrative things related to file permissions.

## Getting started

The data files required for this workshop can be found on our website. You need to download this and place it in your home directory (`/home/username`), before you start this exercise. Please use the commands below to get started (you will learn what these command does later in the exercise). Open the terminal and enter these commands (commands are case sensitive) and each command should be entered in a single line followed by &#8629; (`Enter`) key.
```
cd ~
curl -O http://www.public.iastate.edu/~arnstrm/WORKSHOP_FILES.tar.gz
```
Once your cursor (command prompt) comes back to the original position, type
```
tar xf WORKSHOP_FILES.tar.gz
ls
```
You should see `WORKSHOP_FILES` listed there.

PS: all materials, including the slides, handout and instructions to set up your computer should be in the folder you downloaded.

## Navigation
This section will introduce you to some basic file/directory navigation and manipulation techniques.

### To know the present location
```
pwd
```
Returns you the present working directory (print working directory)

You should see the output something like `/home/username` This means, you are now working in the username directory, which is located in home directory. The directory that you will be in after logging in is your home directory. You can also avoid writing the full path by using ~ in front of your username or simply `~`

`~` or `~username` 	same as	`/home/username`

Present directory is represented as `.` (dot) and parent directory is represented as `..` (dot dot).


### Changing directories
To jump from one directory to another we use the cd (change directory) command.
```
cd ..
```
Changes your present location to the parent directory
```
cd DIRECTORY
```
This changes your location back to your DIRECTORY.

***Task 1.1: Change your directory to the WORKSHOP_FILES directory present in your home directory.***

**TIP**: You can type in first few letters of the directory name and then press `Tab` to auto complete rest of the name (especially useful when the file/directory name is long). This only works when there are unique matches for the starting letters you have typed. If there is more than one matching files/directories, pressing `Tab` twice will list all the matching names. You can also recall your previous commands by pressing &#8593; or &#8595; arrow keys or browse all your previously used commands by typing `history` on your terminal (typically, last 500 commands will be saved in this file).

## Directories and files

### Making directories

To create a directory, `mkdir` (`m`a`k`e `dir`ectory) can be used.

```
mkdir DIRECTORY
```
Unlike PC/Mac folders, here you can’t have space in your directory name (but some special characters are okay). You can also specify the path where you want to create your new folder.

***Task 1.2: Make a new directory named `FirstDirectory` within the `WORKSHOP_FILES` directory. Then change your directory to the `FirstDirectory`.***

```
mkdir FirstDirectory
```

### Copying directories

To copy a file, `cp` (`c`o`p`y) command is used. When using this command you have to provide both source file and destination file.
```
cp SOURCE DESTINATION
```

You can also specify the absolute path of the source and/or destination file. To know more about any command you can use man command, which opens the manual of the command you ask (referred as `man page`).

```
man cp
```
This opens the manual for the `cp` command. Take a look at the manual of `cp` command (use arrow keys to move top or bottom of the page). `OPTIONS` are arguments that can be used to accomplish more from the same command (and are not required for regular operation). Eg., by using option `–i` with the regular `cp` command, you can always make sure that you are not overwriting the existing file while copying (if the target already contains the same file). The syntax for using the options will also be provided in the manual. **To `exit`, press `q`***.

*Looking at the man page for `cp` command, what options can be used to copy a directory (including all files within it)?*

*How else you can get help on cp command (other than ‘man’)?*


***Task 1.3: Now change your directory back to the home directory. Create a copy of `WORKSHOP_FILES` and name it as `BACKUP_WORKSHOP`). This will serve as a backup copy of all files that are required for the workshop (in case you accidentally modify the contents while working).***
```
cp -r WORKSHOP_FILES BACKUP_WORKSHOP
```
