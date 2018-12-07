---
title: "Introduction to Data Acquisition"
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Getting Data from iPlant/Cyverse

Quick tutorial to download data from iPlant datastore to Lightning3/Condo using iRODS. iRODS is currently installed on Lightning3/Condo, to start using it, just load it using modules.
```
module use /data004/software/GIF/modules
module load iRODS
```
iRODS provides Unix command-line utilities to interact with the iPlant Data Store. Many commands are similar to Unix (by adding <code>i</code> to the common UNIX commands, you get iRODS commands or icommands).
### 1. Login to your iPlant account ###

Sometimes also referred as iRODS account, this is required to access private files from the iPlant data store. This is done by logging in via `iinit` command.
```
$ iinit
One or more fields in your iRODS environment file (.irodsEnv) are
missing; please enter them.
Enter the host name (DNS) of the server to connect to:data.iplantcollaborative.org
Enter the port number:1247
Enter your irods user name:aseetharam
Enter your irods zone:iplant
Those values will be added to your environment file (for use by
other i-commands) if the login succeeds.
Enter your current iRODS password:
#password won't be displayed
```
Now you are logged in. Although your working directory doesn't change, all `iCommands` can now access the iPlant datastore.
You can `logout` anytime by typing `iexit` command.

Alternatively, you can save this information in `~/.irods/irodsEnv` file as follows:
```
irodsHost data.iplantcollaborative.org
irodsPort 1247
irodsUserName <username>
irodsZone iplant
irodsPassword <your password>
```
This will let you connect automatically each time you want to login by just typing in `iinit` command. You can also exclude the last line (password), so that it will prompt only password for logging in.

### 2. Download files ###

Navigate to the desired location/directory and use `iget` command to download the desired file. You can also use any of the standard UNIX command by using `i` as prefix (eg., `icp` for copy, `imv` to move, `ipwd` to print working directory etc.,)
```
iget path/to/filename
```
 If there are large number of files to be copied or you just want to sync 2 directories (local and iplant), you can also use `rsync` command. `rsync`
synchronize the data between a local copy (local file system) and the copy stored in iRODS or between two iRODS copies. You can use `rsync` to either:
- synchronization of data from the client's local file system to iRODS
- from iRODS to the local file system
- from one iRODS path to another iRODS path

To synchronize (recursively) the data from the local directory foo1 to the iRODS collection foo2
```
irsync -r foo1 i:foo2
```
To synchronize (recursively) the data from the iRODS collection foo1 to the local directory foo2
```
irsync -r i:foo1 foo2
```
To synchronize (recursively) the data from the iRODS collection foo1 to another iRODS collection foo2
```
irsync -r i:foo1 i:foo2
```

Finally, if you want to transfer a single file from local host to remote (iplant), you can do it using `iput` command.

```
iput filename
```

### References ###
For more information, follow these links:

* [iRODS](https://www.irods.org/index.php/icommands|irods.org)
* [iPlant](https://pods.iplantcollaborative.org/wiki/display/start/Using+icommands|pods.iplantcollaborative.org)


---

[Next](sra.md){: .btn  .btn--primary}
[Previous](file-transfer-using-globus-connect-personal-gcp.md){: .btn  .btn--primary}
[Table of contents](../dAc_introduction.md){: .btn  .btn--primary}
