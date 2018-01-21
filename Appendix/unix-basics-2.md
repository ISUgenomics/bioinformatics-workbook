# Unix Basics 2
This exercise will provide you details about some administrative commands with examples. Here, you can learn how to change permissions for files and folders to modify its accessibility.
### Changing permissions

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
#### 1. Adding permissions

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


#### Removing permissions:

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
```
OPTIONS include `-R` recursively (the permissions are applied to all the files, directories present inside the directory)


***Task 1.12: Check the permissions for the files located in the tutorials directory.***
```
ls -l
```
*What permissions does the group have on these files? Which group does your account belong to?*
