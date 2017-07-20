# Downloading files using wget

Normally, for downloading files we use wget/curl and paste the link (to the file) to download
```
wget http://link.edu/filename
```
But some cases, you may want to download a entire directory of files and downloading directory using wget is not straighforward. 

## Using Wget ##

There are 2 options. You can either specify a regular expression for a file or put a regular expression in the URL itself.
First option is useful, when there are large number of files in a directory, but you want to get only specific format of files (eg., fasta)
```
wget -r --no-parent -A 'bar.*.tar.gz' http://url/dir/
```

The second option is useful if you have numerous files that have the same name, but are in different directory

```
wget -r --no-parent accept-regex=/pub/current_fasta/*/dna/*dna.toplevel.fa.gz ftp://ftp.ensembl.org
```

The files won't be overwritten (as they all have same names), instead they are saved as-is maintaining the directory structure.

Some times, if you have a series of files to download (and are numbered accordingly), you can use UNIX <blockcode> brace expansion</blockcode> 

```
wget http://localhost/file_{1..5}.txt
# this will download file_1.txt, file_2.txt, file_3.txt, file_4.txt and file_5.txt
```

To archive the entire website (yes, every single file of that domain), you can use the mirror option.

```
wget --mirror -p --convert-links -P ./LOCAL-DIR WEBSITE-URL
```

###  Other options to consider  ### 

| Option | What it does | Use case |
| --- | --- | --- |
| `-limit-rate=20k` | Limits Speed to 20KiB/s | Limit the data rate to avoid impacting other users' accessing the server. |
| `-spider` | Check if File Exists | For if you don't want to save a file but just want to know if it still exists. |
| `-w` | Wait Seconds | After this flag, add a number of seconds to wait between each request - again, to not overload a server. |
| `-user=` | Set Username | wget will attempt to login using the username provided. |
| `-password=` | Use Password | wget will use this password with your username to authenticate. |
| `-ftp-user=` or `-ftp-password=` | FTP Credentials | Just like the previous settings, wget can login to an FTP server to retrieve files. |

##Citations  ##

- [Stack-exchange thread](http://unix.stackexchange.com/questions/117988/wget-with-wildcards-in-http-downloads)
- [Blog](http://blog.alastair.pro/2012/10/21/wget-regex-filter-by-file-type/)
- [Forum](http://www.linuxquestions.org/questions/linux-newbie-8/wget-with-regular-expressions-846368/)