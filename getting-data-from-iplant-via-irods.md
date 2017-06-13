<p> Quick tutorial to download data from iPlant datastore to Lightning3/Condo using iRODS. iRODS is currently installed on Lightning3/Condo, to start using it, just load it using modules.</p>
<pre>
module use /data004/software/GIF/modules
module load iRODS
</pre>
<p>iRODS provides Unix command-line utilities to interact with the iPlant Data Store. Many commands are similar to Unix (by adding <code>i</code> to the common UNIX commands, you get iRODS commands or icommands).</p>
<h3> 1. Login to your iPlant account </h3>

<p>Sometimes also referred as iRODS account, this is required to access private files from the iPlant data store. This is done by logging in via <blockcode>iinit</blockcode> command.</p>
<pre>
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
</pre>
<p>Now you are logged in. Although your working directory doesn't change, all <blockcode>iCommands</blockcode> can now access the iPlant datastore.
You can <blockcode>logout</blockcode> anytime by typing <blockcode>iexit</blockcode> command.</p>

<p>Alternatively, you can save this information in <blockcode>~/.irods/irodsEnv</blockcode> file as follows:</p>
<pre>
irodsHost data.iplantcollaborative.org
irodsPort 1247
irodsUserName <username>
irodsZone iplant
irodsPassword <your password>
</pre>
<p>This will let you connect automatically each time you want to login by just typing in <blockcode>iinit</blockcode> command. You can also exclude the last line (password), so that it will prompt only password for logging in.</p>

<h3> 2. Download files </h3>

<p>Navigate to the desired location/directory and use <blockcode>iget</blockcode> command to download the desired file. You can also use any of the standard UNIX command by using <blockcode>i</blockcode> as prefix (eg., <blockcode>icp</blockcode> for copy, <blockcode>imv</blockcode> to move, <blockcode>ipwd</blockcode> to print working directory etc.,)</p>
<pre>
iget path/to/filename
</pre>
<p> If there are large number of files to be copied or you just want to sync 2 directories (local and iplant), you can also use <blockcode>rsync</blockcode> command. <blockcode>rsync</blockcode> 
synchronize the data between a local copy (local file system) and the copy stored in iRODS or between two iRODS copies. You can use <blockcode>rsync</blockcode> to either:</p>
  - synchronization of data from the client's local file system to iRODS
  - from iRODS to the local file system
  - from one iRODS path to another iRODS path

<p>To synchronize (recursively) the data from the local directory foo1 to the iRODS collection foo2</p>
<pre>
irsync -r foo1 i:foo2
</pre>
<p>To synchronize (recursively) the data from the iRODS collection foo1 to the local directory foo2</p>
<pre>
irsync -r i:foo1 foo2
</pre>
<p>To synchronize (recursively) the data from the iRODS collection foo1 to another iRODS collection foo2</p>
<pre>
irsync -r i:foo1 i:foo2
</pre>

<p>Finally, if you want to transfer a single file from local host to remote (iplant), you can do it using <blockcode>iput</blockcode> command.</p>
<pre>
iput filename
</pre>

<h3> References </h3>
For more information, follow these links:

   - [iRODS offical documentation ](https://www.irods.org/index.php/icommands|irods.org)
   - [iPlant basic usage guide](https://pods.iplantcollaborative.org/wiki/display/start/Using+icommands|pods.iplantcollaborative.org)