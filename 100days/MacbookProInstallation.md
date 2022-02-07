---
title: "Macbook Pro Installation"
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

{% include toc %}


## Introduction
This guide is meant for the installation and setup of an M1 MacBook Pro for bioinformatics data analysis. I assume that this is a brand new MacBook and that no software has been installed beyond what is installed with IOS Monterey. A new MacBook is not required for the installation of this software but a clean installation of `conda` and `brew` is recommended.

Setting up your new Mac computing machine is an investment of time that you will quickly appreciate in your daily work. With the help of this tutorial, the process will be straightforward and successful.

## Various Methods of Installation

Depending on how the software is distributed and how popular it is different installation methods are suitable. This guidebook will teach you how to install applications to MacBook Pro that are available in the Apple store or downloaded directly from the distributor's website from the Internet. You will also learn how to install various packages using commands in the terminal. For those who are part of the academic community, in particular Iowa State University, we will also show you how to take advantage of an academic license for software that usually requires a paid license.

### Installation via App Store

The **App Store** is an in-built store platform available on your MacBook Pro that allows you to browse and install the software. The App Store shortcut &nbsp;<img src="https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcS61Gy0FdwGsIx-u4Zu50BPt9ZxbFmsAbRNFApRrs1-zxoQnpAj6VnvXJ3iY3trykjn3Uc&usqp=CAU" alt="Mac App Store" height="18" width="18">&nbsp; is usually visible on your Dock (by default, it is a horizontal menu with icons of pre-installed apps). You can also easily use *Finder* &nbsp;<img src="https://images.macrumors.com/t/5BiCx6nBBb0fGUFWfLHjqaD1zFk=/1200x1200/smart/article-new/2018/02/macos-finder-icon.jpg" alt="Mac App Store" height="18" width="18">&nbsp; or *Launchpad* &nbsp;<img src="https://uploads-ssl.webflow.com/5f7081c044fb7b3321ac260e/5fedaca4acad015c2de7b6a1_30_launchpad.png" alt="Mac App Store" height="22" width="22">&nbsp; to find the App Store on your system. Once you open the app, a search box is located in the upper left corner. After entering your query key, the matching applications will appear. Many of them are free. Then, click the <span style="background-color:#f2f2f2; color:blue;">&nbsp;get&nbsp;</span> button and then the <span style="background-color:green; color:white;">&nbsp;install&nbsp;</span> button. The newly installed program will appear in the Applications folder, easily accessible via the menu on the left side of the *Finder*. The latter is a graphical interface for browsing the disk contents.

### Installation with Dialog Box

Some programs are available for download on the distributor's website as an installation package. This is usually a single file with the extension *.dmg* (or *.pkg*, *.zip.download*, etc.), which can be run by double-clicking after downloading. Then a graphical interface appears. The user is asked to accept the license terms through a dialog box, and the installation process runs almost automatically. Often at the end, the user is asked to manually drag the new program into the Applications folder, which is a good thing to do for order and easy access.

### Installation using Terminal

The vast majority of packages and development modules are available for installation using the command line in the terminal. Stable versions of various programs can be automatically downloaded and installed using helper programs such as `brew`.
Once you get familiar with a coding terminal, you can easily install development packages. Many research groups release their software on code hosting platforms such as [GitHub](https://github.com) or [BitBucket](https://bitbucket.org/). Such repositories can be downloaded and versioned using [Git](https://git-scm.com).


## Basic Software for Easy Daily Work

### Office Productivity Software

Written communication is required in almost every informatics position. The majority of individuals that we will be communicating results will be using either Microsoft Office Word or Apple's Pages.

* Apple Suite (Keynote, Numbers, Pages) should be pre-installed on your MacBook Pro.
* Microsoft Office Suite is available at https://softwarekeep.com/office-for-mac. The Iowa State community can download the entire toolkit with a free academic license through Okta Application Dashboard, available at https://web.iastate.edu/signons. Once the package is downloaded, the automatic installer will guide you through the process. New tools will be available in the Applications on your Mac. Along with the MS package you will be equipped with office tools:
  * **Word** - a word processing software, which offers DOC/DOCX formats,
  * **Excel** - a spreadsheet creator, which features calculation or computation capabilities, graphing tools, and pivot tables,
  * **Powerpoint** - a presentation and slideshow maker,
  * **Teams** - business communication platform, offering workspace chat and videoconferencing, file storage, and application integration,
  * **OneNote** - a note-taking program for free-form information gathering and multi-user collaboration.


### Adobe for Graphic and Design

Adobe Inc. delivers comprehensive solutions for digital media, including Video, Design, Photography, and the Web. The applications can be installed from the Creative Cloud platform. Iowa State has a license for these apps in Okta Application Dashboard, available at https://web.iastate.edu/signons. Otherwise, it is ~$30/month for students/faculty.

Along with the Adobe package you will be equipped with professional tools dedicated for:

* **Photoshop** - editing Photos
* **Illustrator** - editing figures
* **Acrobat DC** - PDF visualization and signing documents
* **InDesign** - multipage graphic/text

`TIPS:` I use the first three more than InDesign.


### Screen Capture & Video Editors

Screenshots can be a quick and effective way to capture important information. For the scientist, it also allows for the preparation of educational materials and quality tutorials that enhanced with graphic visualizations facilitate the assimilation of knowledge by the recipients.
<br>If your intent is to generate tutorials or create video content of your results then you will want some kind of screen capture and video editing software. For very generic screen capture you can use the built-in screen capture or download [QuickTime](https://support.apple.com/guide/quicktime-player/welcome/mac). But if you are planning on generating more professional videos more frequently, I would recommend the following software.
* [Camtasia](https://www.techsmith.com/download/camtasia/), all-in-one screen recorder and video editor - available in the App Store for a fee
* [Apple Pro Apps bundle for education](https://www.apple.com/us-edu/shop/product/BMGE2Z/A/pro-apps-bundle-for-education):
  * Final Cut pro
  * Compressor
  * Motion
  * Logic Pro
* **iScreen Shoter**, screen capture & text scanner - freely available in the App Store
* and many more available in App Store for free (try searching for *'screen capture', 'video record', etc.*)


### Team communication

Virtual meetings, video conferences, and remote collaboration have become our new daily routine. To keep the work efficient and the cooperation strong, we need tools that facilitate live interaction regardless of the distance.
<br>Everyone has their preferences when it comes to remote communication with their team. I have been using Slack with my group since 2016. With a group of experienced bioinformaticians, this form of communication is fantastic as users can very quickly share their knowledge and answer each other's questions. However, if you are the only expert in the group then it can be counterproductive as it becomes a direct conduit to ask you questions. I encourage the members of my group to try to figure out a problem on their own for an hour or two and if they still can't figure it out then to send a slack message to the group to see if anyone else has a solution.
<br>There is a wide range of free tools available for remote video communication and chatting.
* [TEAMS](https://www.microsoft.com/en-us/microsoft-teams/group-chat-software) - included in the MS Office Suite
* [Slack](https://slack.com/downloads/mac) - freely available in the App Store
<br>`TIPS:` The easiest way to reload your teams is to sign in on the slack website and then use your email to sign in for each team you want to be signed in on your desktop application.
* [Zoom](https://zoom.us/support/download?os=mac) - freely available in the App Store
<br>`TIPS:` Don't forget to copy over your zoom backgrounds
* [Webex](https://www.webex.com) - freely available in the App Store


### Password Manager

A password manager is strongly recommended. We have passwords now for everything from multiple emails (work, home) to University-related tasks (e.g., *Workday*) to our favorite social media (*Twitter, Facebook, LinkedIn*) to our favorite way to relax (*Apple TV, Netflix, Amazon Prime, etc.*). Trying to remember all of these passwords is an impossible task unless you start reusing passwords which is a bad idea with the number of companies getting information hacked.
<br>A password manager is a computer program that allows users to store, generate, and manage their passwords for local applications and online services. A password manager keeps everything organized, easy to find, and can generate long random passwords that are essentially impossible to break. With this solution, your passwords are encrypted and digitally secure. Consider the following options for your new MacBook Pro:
* **[Dashlane](https://www.dashlane.com)** - freely available in the App Store
* **[OnePass](https://1password.com/)** - freely available in the App Store

`TIPS:` I personally use Dashlane, but I have heard great things about OnePass as well.


### Web Browser Extensions

For Apple devices, Safari is the default built-in browser &nbsp;<img src="https://purepng.com/public/uploads/large/purepng.com-safari-iconsymbolsiconsapple-iosiosios-8-iconsios-8-7215225961106timx.png" alt="Mac App Store" height="22"  width="22">&nbsp;. You can also install other browsers yourself, e.g. here you can find detailed [installation instructions](https://www.techsolutions.support.com/how-to/how-to-download-and-install-google-chrome-on-a-mac-12424) for Google Chrome.
<br>`TIPS:` If you are a web-developer, it is a good idea to have several different browsers to pre-test your implementations.

Have too many tabs open at once and don't want to close them just in case you decide you will come back sometime later? Well, I have that problem too, and it slows down your computer. Instead, install the browser extension and have it save the entire window of tabs in a list that you can resummon with a click.
For Safari it is useful to add a few extensions:

 * **OneTab**, allows you to save your session for later or reduce all tabs into one list - freely available in the App Store


### Productivity Software

 I have found the following programs really useful to streamline my workflows and better interact with my laptop and mouse.

 * [Better Touch Tool](https://folivora.ai) - available as *BetterSnapTool* in the App Store for a small fee
   * software to customize different input devices like your mouse and trackpad and window snapping.
 * [Alfred](https://www.alfredapp.com) - freely available in the App Store
   *  hotkeys, keywords, text expansion and more


## Install Basic Developer Tools

### Install Xcode
Xcode is Apple's integrated development environment for macOS, used to develop software. In order to get the most out of your MacBook, this developer's toolkit is required. Installing XCode at a minimum also installs git for *GitHub* commands. There are probably ways where you don't even have to install XCode developer tools and still get *GitHub* and the other functionality, but XCode installation seems to work well. You will mostly use the command line tools unless you want to learn swift and ios app development.

* The download takes a long time
* The unzipping takes a long time
* Move the Xcode App to your applications

Install the XCode Command line tools:

```
xcode-select --install
```

### Install Rosetta 2
For those applications that are still compiled for Intel processors, Apple has Rosetta 2, which is a translation environment for your favorite intel compiled programs.

Run the following command on the terminal:

```
/usr/sbin/softwareupdate --install-rosetta --agree-to-license
```

Then, restart your computer.


### Install Homebrew

Homebrew is a package manager that makes installing many useful packages really easy. In short, it installs the software you need that Apple didn’t. The tool installs packages into their own directory and then symlinks their files into `path:/usr/local` on macOS.

###### Download HomeBrew:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
```

Be sure to change *'YourNAME'* to your local username.
```
echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> /Users/YourNAME/.zprofile
   eval "$(/opt/homebrew/bin/brew shellenv)"

```

###### Create a brew file:

I update this brew file as I add more programs with brew. That way, I will always know what I have installed and can quickly install it on a new machine.

Place the following text into a file called `brewfile`.

```
tap "homebrew/bundle"
tap "homebrew/core"
tap "homebrew/cask"
brew "zsh-completions"
brew 'bash-completion'
brew "gnu-sed"
brew "gnu-tar"
brew "gawk"
brew "gnu-which"
brew "gzip"
brew "unzip"
brew "coreutils"  #Adds a few extra commands typically found in Unix systems
brew "curl" #Updated curl
brew "dos2unix" #For those pesky dos line endings
brew "findutils" #find xargs locate updatedb
#brew "git" #if you installed XCODE then this is already installed
brew "go" #programming language
brew "jq" #https://stedolan.github.io/jq/ #like sed for JSON
brew "ruby", link: true #programming language
brew "rust" #programming language
brew "tree" #see your file structure
brew "wget" #like curl but better
cask "docker" #Not singularity but can be useful.
cask "iterm2" #An advanced Terminal.
cask "xquartz" #Required for some plotting programs like R
cask "figtree" #Visualize phylogenic trees
brew "tcl-tk" #Needed for some programs
tap "brewsci/bio/"
brew "brewsci/bio/pymol" #Visualize protein structures
brew "igv"  #Genome browser
cask "jbrowse" #A Better Genome browser
brew "htop"  #A different type of top for your mac
brew "pygments" #color syntax

brew "fastqc"
```

Execute the following command in the same folder as the `brew file` defined above and it will install all of the programs.

```
brew bundle
```

### Install Atom editor
Having a great markdown editor to go along with your *GitHub* repo is a must for documentation of bioinformatic projects.
<br>**Atom** is a free and open-source text and source code editor available for download at https://atom.io. This editor should be at the core of every developer’s toolbox. It allows for cross-platform editing, smart autocompletion for multiple programming languages, find-preview-replace actions, browsing file system, editing code in multiple panes, managing built-in packages, and versioning code using Git.

Start the installation process by double-clicking on the file you downloaded to your Downloads folder.
<br>`TIPS:` Atom did not open for me the first time when I double-clicked on it. You may have to right-click on the application and select *'Open'*.

Here are some recommended packages, which can be easily incorporated into the editor's functionality by installing them using the built-in package manager.

* language-swift-89
* language-r
* markdownn-folding
* markdown-pdf
* minimap
* wordcount
* drag-relative-path
* markdown-scroll-sync
* autocomplete-python
* autocomplete-swift
* autocomplete-R


### Terminal setup

###### *Install iterm2*

It is a powerful terminal. I haven't utilized its features fully. One option I felt was missing was 'skip by the word' on the command line, but apparently, that is a really-easy fix since iterm2 is fully configurable.  

[How to skip by word in iterm2](https://coderwall.com/p/h6yfda/use-and-to-jump-forwards-backwards-words-in-iterm-2-on-os-x)
  * esc + f  
  * esc + b

###### *Add Oh-My-ZSH to make the terminal more useful*

```
sh -c "$(curl -fsSL https://raw.githubusercontent.com/robbyrussell/oh-my-zsh/master/tools/install.sh)"
```

This generates a new .zshrc file.  Change the line with plugins to the following.  Include what is useful to you.  I don't use sublime or vscode but they are quite popular.

```
plugins=(git z github history osx pip pyenv pylint python sublime vscode)
```

###### *Setting up .zshrc / .bashrc / .bash_profile*

```
# Aliases
pcat='pygmentize -f terminal256 -O style=native -g'
```


###### *Giving Terminal Full Disk Access*

If you start using Terminal, Mac will now ask you every time if you change directories permission to change into that directory.  To give Terminal access to all folders:

* Settings
  * Security & Privacy
    * Full Disk Access


###### *Install midnight commander file menager*

Midnight Commander (**mc**) is a free cross-platform visual file manager that can be used as a primary file manager in a terminal session. By default, it consists of two panels, allowing simultaneous viewing from a terminal of two locations in the file system. That can also be a remote location (e.g., sshfs mounted HPC cluster), which allows you to view graphical files with your favorite programs without downloading the files to your local machine.

Install midnight-commander with brew:
```
brew install midnight-commander
```
And run it calling `mc` command in the terminal. Use [mc cheat sheet](https://gist.github.com/sgergely/3793166) for Mac to familiarize yourself with keyboard shortcuts.


### Git setup

You will want to copy over the following files from your old machine to your new machine to make github and Atom work again with pushes and pulls.  Github does not allow username/password pushes/fetches anymore.

* .gitconfig
* .git-credentials

I also had to create a new ssh key in my github account for my new laptop. If you need a refresher on how to do this see this [GitHub tutorial](https://bioinformaticsworkbook.org/Appendix/github/introgithub.html#gsc.tab=0).



### Install VPN client

Virtual Private Network (VPN) is a technology that encrypts your internet traffic improve your online privacy. This means the data you transfer between your MacBook Pro and the HPC infrastructure is protected. Many computing networks, including SCINet and ISU HPC, require a secure connection via VPN to access resources and schedule your computations.

If you need a VPN to access your HPC you will need to download the software suggested by the administrator of the resource you wish to access. For us at Iowa State University, we use the [Cisco AnyConnect client](https://www.cisco.com/). The details of the installation and your first connection can be found at https://it.engineering.iastate.edu/how-to/install-and-connect-to-vpn-pc/.


### File Transfer software

File transfer is the translocation or transport of data, usually stored in the file, through a communication channel from one computer system to another. Every time you connect via VPN to a computing infrastructure and want to put your input there, you will be transferring files. If your data is large, it's worth equipping yourself with a tool that will make file transfer easier, faster, and more secure. Consider the following options for your new MacBook Pro:
* [filezilla](https://filezilla-project.org/download.php)
  * Copy over any site information you need before formatting your old computer.
* [Globus](https://www.globus.org/data-transfer)
  * Move, sync, and share large amounts of data between GridFTP server and a user's computer.
* **sshfs**
  * Mount a remote filesystem from another server locally on your machine using SFTP. This is very useful for viewing files held on the cluster with your local graphics programs (e.g. charts, text files both .docx and .pdf).

###### Install SSHFS

Download via terminal both `macfuse` and and `sshfs` from osxfuse GitHub repository.

```
wget https://github.com/osxfuse/osxfuse/releases/download/macfuse-4.2.4/macfuse-4.2.4.dmg
wget https://github.com/osxfuse/sshfs/releases/download/osxfuse-sshfs-2.5.0/sshfs-2.5.0.pkg
```
Then, find the files in your Downloads folder and run the executable files in the same order.


## Install Developer Libraries

### Install R

Install R from the package from R cran site.

* [Install R](https://mirror.las.iastate.edu/CRAN/)
  * https://mirror.las.iastate.edu/CRAN/
    * https://mirror.las.iastate.edu/CRAN/bin/macosx/big-sur-arm64/base/R-4.1.2-arm64.pkg

`TIPS:` I got the following error when I tried to install from DMG so don't do that.
**Error** "R can't be opened because Apple cannot check it for malicious software"
  * [Monterey Arm version 4.2](https://mac.r-project.org/monterey/R-devel/R-GUI-8008-4.2-monterey-arm64-Release.dmg)

### Install R Studio
Up your R desktop game by installing R studio.

* https://www.rstudio.com/products/rstudio/download/#download


#### *Check to see R is working*

I use plot to quickly test the installation.

```
plot(1:10)
```

A window should pop up with a scatter plot along the diagonal.


### Install R Manager [Renv](https://cran.r-project.org/web/packages/renv/)

I haven't explored this enough yet, but it would make sense to have a good package manager for all the potential R packages you might install. Since some can conflict, it can also be good to have separate conda environments for R packages.

* https://cran.r-project.org/web/packages/renv/
* https://6chaoran.wordpress.com/2020/07/20/introduction-of-renv-package/


pacman?
R environment manager

mamba conda

update



### Install Python Modules


```
pip install multiqc
```


### Install Conda

Conda is an open-source and cross-platform environment management system. Conda quickly installs, runs, and updates packages and their dependencies. It helps build the virtual environment required for a specific project, including novel software development.
<br>Conda is yet another package manager that is very popular in the bioinformatics community. Almost every software you want to install can be installed using Conda by creating a `conda` environment. The new MacBook Pros with the M1 Arm chips does make installations a little more challenging as not all software has been formatted to run natively on the ARM architecture. Fortunately, Apple's translation environment software, Rosetta 2, can be used to install and run anything we want that was meant for Intel chips.

Use the [conda cheat sheet](https://kapeli.com/cheat_sheets/Conda.docset/Contents/Resources/Documents/index) to get an overview of the available functionalities.


###### Install Miniforge3 for both ARM and Intel chips.


* https://github.com/conda-forge/miniforge
  * [Miniforge3](https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh)
  * [Miniforge3_x86](https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh)

I installed the x86 version in a folder with `_x86` at the end of it. This will become important later.

```
                          ~/Software/miniforge3
base                  *   ~/Software/miniforge3_x86
```

This now gives us two base conda environments.  One for installations native to the ARM architecture and native to M1 macs and one for all the other programs that haven't made an ARM version. The x86 programs will run using `rosetta2` translation environment.

###### *How to change between base conda installations*
This website does a really good job explaining that we just need to change the code in the `.zshrc` file: [Changing base conda installs](https://stackoverflow.com/questions/58131555/how-to-change-the-path-of-conda-base).

I placed all of the next code at the very end of this file. Oh-My-ZSH has a lot of other text in this file.  Leave that alone.

* .zshrc

This file of course will look slightly different as you will have placed the miniforge3 folder in a different location then on my laptop.

```bash
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/Users/andrewseverin/GIFNew/Software/miniforge3/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/Users/andrewseverin/GIFNew/Software/miniforge3/etc/profile.d/conda.sh" ]; then
        . "/Users/andrewseverin/GIFNew/Software/miniforge3/etc/profile.d/conda.sh"
    else
        export PATH="/Users/andrewseverin/GIFNew/Software/miniforge3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

# this provides the prompt the conda env name you are in.  OH-My-ZSH doesn't do this automatically unfortunately.
POWERLEVEL9K_LEFT_PROMPT_ELEMENTS=(anaconda ...ENVS)
```

The main point that the website above makes is that in order to change the base installation, all we have to do is change the folder name from `miniforge3` to `miniforge3_x86`.  Doing this every time we want to change between base installations would be a real pain so I modified the script to make this a lot easier.


```bash

if [[ $x86true -eq 1 ]]
then
x86="_x86"
else
x86=""
fi

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/Users/andrewseverin/GIFNew/Software/miniforge3$x86/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/Users/andrewseverin/GIFNew/Software/miniforge3$x86/etc/profile.d/conda.sh" ]; then
        . "/Users/andrewseverin/GIFNew/Software/miniforge3$x86/etc/profile.d/conda.sh"
    else
        export PATH="/Users/andrewseverin/GIFNew/Software/miniforge3$x86/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

# this provides the prompt the conda env name you are in.  OH-My-ZSH doesn't do this automatically unfortunately.
POWERLEVEL9K_LEFT_PROMPT_ELEMENTS=(anaconda ...ENVS)
```

As you can see I added an if statement that changes a variable `$x86` which I placed at the the end of the folder name to modify the conda location. With this modification, we can create two new functions that will permit us to very quickly change between base conda installations.


* condaArm

```
x86true=0
source ~/.zshrc
conda info --envs
```

* condaX86
```
x86true=1
source ~/.zshrc
conda info --envs
```

To change between them all you have to do is source the file

```
source condaX86
```

or

```
source condaArm
```



## Install Custom Packages


### Qiime2

QIIME2 is a next-generation microbiome bioinformatics platform. If you do not run metagenomic analyses yet, then you can skip this setup.
<br>We do a lot of metagenomic analyses, and Qiime can be one of the more challenging software to install. In this tutorial, I will use it as an example of how to install software on your local machine using the conda environment.

* [Qiime install directions](https://docs.qiime2.org/2021.11/install/native/#miniconda)

```
source ~/condaX86

conda update conda
wget https://data.qiime2.org/distro/core/qiime2-2020.11-py36-osx-conda.yml
```

Here I am installing mamba for X86. As an install of mamba will initially run only for the ARM architecture, I had to create a specific environment to install the X86 version of mamba. So, I change into it to use it to install X86 programs.

The only downside to this is that I have to use the full path in order to source active the conda environment.
```
conda create --name mambaX86
conda install -n mambaX86 -c conda-forge mamba
mamba env create -n qiime2-2021.11 --file qiime2-2021.11-py38-osx-conda.yml
```

To run Qiime just activate the conda environment.

```
# you have to use your full path in this case.
conda activate /Users/andrewseverin/GIFNew/Software/miniforge3_x86/envs/mambaX86/envs/qiime2-2021.11

```
