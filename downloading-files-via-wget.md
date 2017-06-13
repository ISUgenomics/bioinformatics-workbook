<p>Normally, for downloading files we use wget/curl and paste the link (to the file) to download</p>
<pre>
wget http://link.edu/filename
</pre>
<p>But you can also download entire file in the directory that matches your regular expression using the examples below</p>

<h2> Using Wget</h2>

<p>There are 2 options. You can either specify a regular expression for a file or put a regular expression in the URL itself.
First option is useful, when there are large number of files in a directory, but you want to get only specific format of files (eg., fasta)</p>
<pre>
wget -r --no-parent -A 'bar.*.tar.gz' http://url/dir/ 
</pre>

<p>The second option is useful if you have numerous files that have the same name, but are in different directory</p>

<pre>
wget -r --no-parent accept-regex=/pub/current_fasta/*/dna/*dna.toplevel.fa.gz ftp://ftp.ensembl.org
</pre>

<p>The files won't be overwritten (as they all have same names), instead they are saved as-is maintaining the directory structure.

Some times, if you have a series of files to download (and are numbered accordingly), you can use UNIX <blockcode> brace expansion</blockcode> </p>

<pre>
wget http://localhost/file_{1..5}.txt
# this will download file_1.txt, file_2.txt, file_3.txt, file_4.txt and file_5.txt
</pre>

<p>To archive the entire website (yes, every single file of that domain), you can use the mirror option.</p>

<pre>
wget --mirror -p --convert-links -P ./LOCAL-DIR WEBSITE-URL
</pre>

<h3> Other options to consider </h3>

<table>
<thead><tr><th>Option</th><th>What it does</th><th>Use case</th></tr></thead><tbody>
 <tr><td><blockcode>-limit-rate=20k</blockcode></td><td>Limits Speed to 20KiB/s</td><td>Limit the data rate to avoid impacting other users' accessing the server.</td></tr>
 <tr><td><blockcode>-spider</blockcode></td><td>Check if File Exists</td><td>For if you don't want to save a file but just want to know if it still exists.</td></tr>
 <tr><td><blockcode>w</blockcode></td><td>Wait Seconds</td><td>After this flag, add a number of seconds to wait between each request - again, to not overload a server.</td></tr>
 <tr><td><blockcode>-user=</blockcode></td><td>Set Username</td><td>wget will attempt to login using the username provided.</td></tr>
 <tr><td><blockcode>-password=</blockcode></td><td>Use Password</td><td>wget will use this password with your username to authenticate.</td></tr>
 <tr><td><blockcode>-ftp-[user|password]=</blockcode></td><td>FTP Credentials</td><td>Just like the previous settings, wget can login to an FTP server to retrieve files.</td></tr>
</tbody></table>
<h2> Citations </h2>

  - [Stack-exchange thread](http://unix.stackexchange.com/questions/117988/wget-with-wildcards-in-http-downloads)
  - [Blog](http://blog.alastair.pro/2012/10/21/wget-regex-filter-by-file-type/)
  - [Forum](http://www.linuxquestions.org/questions/linux-newbie-8/wget-with-regular-expressions-846368/)