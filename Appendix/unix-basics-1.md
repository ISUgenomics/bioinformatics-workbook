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
