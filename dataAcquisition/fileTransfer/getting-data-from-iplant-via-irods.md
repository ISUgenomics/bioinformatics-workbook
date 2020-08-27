---
title: "Introduction to Data Acquisition"
layout: single
author: Arun Seetharam
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---



# Getting Data from iPlant/CyVerse

This is the quick tutorial to download data from [CyVerse](https://en.wikipedia.org/wiki/IPlant_Collaborative) (previously known as iPlant) on HPC clusters using [iRODS](https://irods.org/). iRODS may or may not be available as module on your cluster, and if available as a module load the module to use `irods` commands.

```bash
# search
module spider irods # or
module spider icommands
# load
module load <module name>
```

If the module is unavailable, you can either install it yourself following the guidelines [here](https://docs.irods.org/master/getting_started/installation/#irods-setup) or skip installation by using the singularity container (and run these commands within the container)

```bash
#pull container
singularity pull --name irods.sif shub://mjstealey/singularity-irods-icommands
# use container
singularity shell --bind $PWD irods.sif
singularity> iinit
```

### 1. Login to your CyVerse account ###

Sometimes also referred as iRODS account, this is required to access private files from the CyVerse data store. This is done by logging in via `iinit` command. You only need this once per cluster. As it involves data transfer, the best practice recommends using the data transfer nodes on HPC.

```
$ iinit
One or more fields in your iRODS environment file (.irodsEnv) are
missing; please enter them.
Enter the host name (DNS) of the server to connect to: data.iplantcollaborative.org
Enter the port number: 1247
Enter your irods user name: <username>
Enter your irods zone: iplant
Those values will be added to your environment file (for use by
other i-commands) if the login succeeds.
Enter your current iRODS password:
#password won't be displayed
```

Now you are logged in. Although your working directory doesn't change, all `iCommands` can now access the CyVerse datastore.
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

A full list of irods commands can be found [here](https://docs.irods.org/master/icommands/user/)


### References ###
For more information, follow these links:

* [iRODS](https://www.irods.org/index.php/icommands)
* [CyVerse](https://pods.iplantcollaborative.org/wiki/display/start/Using+icommands)


---

[Next](sra.md){: .btn  .btn--primary}
[Previous](file-transfer-using-globus-connect-personal-gcp.md){: .btn  .btn--primary}
[Table of contents](../dAc_introduction.md){: .btn  .btn--primary}
