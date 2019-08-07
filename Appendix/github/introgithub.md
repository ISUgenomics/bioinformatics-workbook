---
title: "Useful Programs and Unix Basics"
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Introduction to GitHub

Github is a code hosting platform for version control and collaboration.  There are three primary use cases for using a version control system like github in science.
  * Sharing of bioinformatics **scripts**
  * Bioinformatic **Program development**
  * **Documentation** of bioinformatic analyses

Probably the most important use-case for new users is documentation.  For those transitioning from a wet-lab, git-hub repos can be thought of as the equivalent to a web-lab notebook, where every command performed in a bioinformatics analysis is recorded with an explanation as to why it was performed, when it was performed (date) and where it was performed (pwd). Github documents can be written using markdown (See [markdown](../Markdown.md) tutorial for an introduction).

---
# How to get a Github account

Signing up for an account is very easy.  Just go to the  [signup webpage](https://github.com/join?source=header-home) and fill out the form.

Most choose the free Unlimited public repositories option and don't set up an organization right away.

After you have an account, and if you are a researcher or educator you can request a free upgrade at [about-github-education-for-educators-and-researchers/](https://help.github.com/articles/about-github-education-for-educators-and-researchers/).

---

# Setting up your Authentication key to streamline remote access

## Authentication with a SSH key

Passwords are not always secure and can be annoying to type.

* SSH keys are much more secure and allow you to log in without typing your password
(or just a different, simpler passphrase).
* When you generate a key, you create two things: a public key and a private key.
* You place the public key on any server and then unlock it by connecting to it with a client that already has the private key.
* When the two match up, the system unlocks without the need for a password.
* SSH keys are also very important for using Git with remote hosts (e.g., GitHub)
---

## Setup Authentication with a SSH key

Log in to the `Remotemachine`

```
ssh <yourID>@remoteMachine
```
or open a terminal on your local machine.

Create the key pair in your home directory:

```
$ ssh-keygen -t rsa
```

Once the `ssh-keygen` command has been issued, you will be asked a few questions:

```
Enter file in which to save the key (/home/yourID/.ssh/id_rsa):
```

You can just hit _enter_ here and it should save it to the file path given.

```
Enter passphrase (empty for no passphrase):
```

 Here is where you decide if you want to password-protect your key. The downside, to having a passphrase, is then having to type it in each time you use the Key Pair.

---
## Create the SSH key

```
$ ssh-keygen -t rsa
Generating public/private rsa key pair.
Enter file in which to save the key (/home/userid/.ssh/id_rsa):
Enter passphrase (empty for no passphrase):
Enter same passphrase again:
Your identification has been saved in /home/userid/.ssh/id_rsa.
Your public key has been saved in /home/userid/.ssh/id_rsa.pub.
The key fingerprint is:
4a:dd:0a:c6:35:4e:3f:ed:27:38:8c:74:44:4d:93:67 userid@Arrow
The key's randomart image is:
+--[ RSA 2048]----+
|          .oo.   |
|         .  o.E  |
|        + .  o   |
|     . = = .     |
|      = S = .    |
|     o + = +     |
|      . o + o .  |
|           . o   |
|                 |
+-----------------+
```

---

## Copy SSH key to GitHub

You now have a file called `id_rsa.pub` in your `.ssh` folder.

```
cat id_rsa.pub
```

```
[userid@hpc-class ~]$ cat .ssh/id_rsa.pub
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQC8QgicqcpFPeyYZhJFW5lBTAdAjHBaYzLwH3l7+lrdmpEKMMMhXMZV5ucxN5WzWU/ERYviYQvQ8NBzkSuHo+SgNJkufF92UqeHIfI/KqgVEGbQn6NGfa5WFBgWZAJAjMzUUrAhJ2qsBez4M1f70os1S2SNcNfoFAJRdWEGE8SX2lpww8+VdCOY6ONw3AYbZbrZtn/ua2hJg+XjYb3T04ggV6TNyV4nnN5r2pRIjJA5OBX1TWcB9HOE4ZIGZoZlk5OYuUJ5rOfuzVLqanWayB3ujuPxW3IUmI6XJt7uDc1N5iVNs2FusjSZmuggWtzCw/pb7EExvNxYMYOxCsewjE0L userid@<remotehost>
```

Copy the entire contents of this file and add it to your GitHub account.

---

# How to Add SSH key to GitHub Repo

Add your new ssh key to your GitHub account by going to [_SSH and GPG keys_](https://github.com/settings/profile/keys) in your profile Settings.  You will need to be logged into GitHub for the above link to work.

![](https://codenvy.com/docs-version/5.0.1/docs/assets/imgs/Clipboard6.jpg)

image source: [https://codenvy.com/docs/user-guide/git-svn/index.html](https://codenvy.com/docs/user-guide/git-svn/index.html)

---
# Version Control Systems

## What do they allow you to do?

* Track changes made to each file
* Revert the entire project or a single file to a previous version
* Review changes made over time
* View who modified the file (and `blame` them for something if necessary)
* Collaborate with others without overwriting their work or risk file corruption, etc.
* Have multiple independent **_branches_** of the same repository and make changes without
effecting others' work.
* And more...

---


# Why to &#x2665; Git?

## Git manages a filesystem as a set of snapshots

Snapshots are called _commits_

![](https://git-scm.com/book/en/v2/book/01-introduction/images/snapshots.png)




(image source:
[https://git-scm.com](https://git-scm.com/book/en/v2/Getting-Started-Git-Basics))

---


# &#x2665; Git &#x2665;

## Almost every interaction with Git happens locally

<div style="text-align:center"><img src="https://git-scm.com/book/en/v2/book/01-introduction/images/areas.png" alt="gource" width="650" /></div>

(image source [https://git-scm.com](https://git-scm.com/book/en/v2/Getting-Started-Git-Basics))


---

# &#x2665; Git &#x2665;

### A remote host adds an additional level to a Git repository

Also, allows for collaboration and back-up.

--

<div style="text-align:center"><img src="https://2.bp.blogspot.com/-w2Z0FGVzygM/WBoAt2kuuCI/AAAAAAAABtU/voWHlUobfl8jB5-5NuSZD0BlN3A9jdMtgCLcB/s1600/basic-remote-workflow.png" alt="Drawing" width="425" /></div>

(image source [https://www.git-tower.com](https://www.git-tower.com/learn/git/ebook/en/command-line/remote-repositories/introduction))

---

Contibutions for this markdown document came from Matt Hufford and was modified with permission from his [BCB546](https://github.com/EEOB-BioData/BCB546X-Fall2017/tree/master/Week_03/) course.  Check out his amazing power point he created at the above link.

---
[Table of contents](../programs.md)
