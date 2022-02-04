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
---

This guide is meant for the installation and setup of an M1 Macbook Pro for bioinformatic data analysis. I assume that this is a brand new macbook and that no software has been installed beyond what is installed with IOS Monterey.  A new macbook is not required for the installation of these software but a clean installation of `conda` and `brew` is recommended.   


## Optional Software

---

#### Word Processing

Written communication is required in almost every informatics position.  The majority of individuals that we will be communicating results will be using either Microsoft Office word or Apple's Pages.

* Microsoft Office Suite (Word, Excel, Powerpoint, Teams, OneNote)
* Apple Suite (Keynote, Numbers, Pages)

#### Password Manager

A password manager is strongly recommended. We have passwords now for everything from multiple emails (work, home) to University related tasks (workday) to our favorite social media (twitter, facebook, linkedIn) to our favorite way to relax (Apple TV, Netflix, Amazon Prime etc).  Trying to remember all of these passwords is an impossible task unless you start reusing passwords which is a bad idea with the amount of companies getting information hacked. A password manager, keeps everything organized, easy to find and can generate long random passwords that are essentially impossible to break.  

* Password Manager (Dashlane, OnePass, etc)

I personally use Dashlane, but I have heard great things about OnePass as well.

#### Screen Capture/Video Editors


If your intent is to generate tutorials or create video content of your results then you will want some kind of screen capture and video editing software.  For very generic screen capture you can use the built in screen capture or download QuickTime.  But if you are planning on generating more professional videos more frequently, I would recommend the following software.

* [Camtasia](https://www.techsmith.com/download/camtasia/)
* [Apple Pro Apps bundle for education](https://www.apple.com/us-edu/shop/product/BMGE2Z/A/pro-apps-bundle-for-education)
  * Final Cut pro
  * Compressor
  * Motion
  * Logic Pro

## Team communication

---

Everyone has their preferences when it comes to remote communication with their team.  I have been using Slack with my group since 2016. With a group of experienced bioinformaticians, this form of communication is fantastic as users can very quickly share their knowledge and answer each other's questions.  However, if you are the only expert in the group then it can be counter productive as it becomes a direct conduit to ask you questions.  I encourage the members of my group to try to figure out a problem on their own for an hour or two and if they still can't figure it out then to send a slack message to the group to see if anyone else has a solution.

* [Slack](https://slack.com/downloads/mac)
  * The easiest way to reload your teams is to sign in on the slack website and then use your email to sign in for each team you want to be signed in on your desktop application.
* [Zoom](https://zoom.us/support/download?os=mac)
  * Don't forget to copy over your zoom backgrounds



## Install Xcode Developer tools

In order to get the most out of your macbook, this developer's toolkit is required.  Installing XCode at a minimum also installs git for github commands. There are probably ways where you don't ever have to install XCode developer tools and still get github and the other functionality but XCode installation seems to work well. You will mostly use the command line tools unless you decide you want to learn swift and ios app development.

* The download takes a long time
* The unzipping takes a long time
* Move the Xcode App to your applications

Install the XCode Command line tools

```
xcode-select --install
```

## Install Rosetta 2
For those applications that are still compiled for Intel processors, Apple has Rosetta 2, which is a translation environment for your favorite intel compiled programs.

Run the following command on the terminal

```
/usr/sbin/softwareupdate --install-rosetta --agree-to-license
```

Restart your computer

## Install Atom editor

Having a great markdown editor to go along with your github repo is a must for documentation of bioinformatic projects.

* https://atom.io

Atom did not open for me the first time when I double clicked on it. You may have to right click on the applicaton and select open.

Here are some recommended packages

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

## VPN client

If you need a vpn to access your HPC you will need to download that.  For us at Iowa State, we use the Cisco Any Connect client

## Adobe

These apps can be installed from the creative cloud app. Iowa State has a license for these in Okta. Otherwise, it is ~$30/month for students/faculty.

* Photoshop
  * For editing Photos
* Illustrator
  * for editing figures
* Acrobat DC
  * for PDF visualization and signing
* InDesign
  * for multipage graphic/text

I use the first three more than InDesign.

# Safari Extensions

Have too many tabs open at once and don't want to close them just in case you decide you will come back sometime later?  Well I have that problem too and it slows down your computer.  Instead install this extension and have it save the entire window of tabs in a list that you can resummon  with a click.

* OneTab (you can download this from App extensions store from Apple)

## Install R

Install R from the package from R cran site.

* [Install R](https://mirror.las.iastate.edu/CRAN/)
  * https://mirror.las.iastate.edu/CRAN/
    * https://mirror.las.iastate.edu/CRAN/bin/macosx/big-sur-arm64/base/R-4.1.2-arm64.pkg

I got the following error when I tried to install from DMG so don't do that.
**Error** "R can't be opened because Apple cannot check it for malicious software"
  <!-- * [Monterey Arm version 4.2 <- don't use](https://mac.r-project.org/monterey/R-devel/R-GUI-8008-4.2-monterey-arm64-Release.dmg) -->

## R Studio
Up your R desktop game by installing R studio.

* https://www.rstudio.com/products/rstudio/download/#download


#### Check to see R is working

I use plot to quickly test the installation.

```
plot(1:10)
```

A window should pop up with a scatter plot along the diagonal.


## R package manager [Renv](https://cran.r-project.org/web/packages/renv/)

I haven't explored this enough yet but it would make sense to have a good package manager for all the potential R packages you might install.  Since some can comflict it can also be good to have separate conda environments for R packages.

* https://cran.r-project.org/web/packages/renv/
* https://6chaoran.wordpress.com/2020/07/20/introduction-of-renv-package/

<!--
pacman?
R environment manager

mamba conda

update -->

## Productivity software

I have found the following programs really useful to streamline my workflows and better interact with my laptop and mouse.

* [Better Touch Tool](https://folivora.ai)
  * software to customize different input devices like your mouse and trackpad and window snapping.
* [Alfred](https://www.alfredapp.com)
  *  hotkeys, keywords, text expansion and more

## File Transfer software

* [filezilla](https://filezilla-project.org/download.php)
  * Copy over any site information you need before formatting your old computer

## Install HomeBrew

Homebrew is a package manager that makes installing many useful packages really easy.  

#### Download HomeBrew

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
```

Be sure to change YourName to your local username.

```

echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> /Users/YourNAME/.zprofile
   eval "$(/opt/homebrew/bin/brew shellenv)"

```

#### Create a brewfile

This brew file I update as I add more programs with brew, that way I will always know what I have installed and can quickly install it on a new machine.

Place the following text into a file called `brewfile`
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

Execute the following command in the same folder as the brewfile defined above and it will install all of the programs.

```
brew bundle
```

## Install python modules

```
pip install multiqc
```



#### Notes about iterm2

It is a powerful terminal and I haven't utilized its features fully. One feature that I felt was missing was being able to skip by word on the command line but apparently that is a really easy fix since iterm2 is fully configurable.  

[How to skip by word in iterm2](https://coderwall.com/p/h6yfda/use-and-to-jump-forwards-backwards-words-in-iterm-2-on-os-x)
  * esc + f  
  * esc + b

## Add Oh-My-ZSH to make the terminal more useful.

```
sh -c "$(curl -fsSL https://raw.githubusercontent.com/robbyrussell/oh-my-zsh/master/tools/install.sh)"
```

This generates a new .zshrc file.  Change the line with plugins to the following.  Include what is useful to you.  I don't use sublime or vscode but they are quite popular.

```
plugins=(git z github history osx pip pyenv pylint python sublime vscode)
```


## Conda

Conda is yet another package manager that is very popular in the bioinformatics community.  Almost every software you want to install can be installed using Conda by creating a conda environment.  The new Macbook Pros with the M1 Arm chips does make installations a little more challenging as not all software has been formatted to run natively on the ARM architecture.  Fortuntely, Apples translation environment software, Rosetta 2, can be used to install and run anything we want that was meant for Intel chips.

Install Miniforge3 for both ARM and Intel chips

* https://github.com/conda-forge/miniforge
  * [Miniforge3](https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh)
  * [Miniforge3_x86](https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh)

I installed the x86 version in a folder with `_x86` at the end of it. This will become important later.

```
                          ~/Software/miniforge3
base                  *   ~/Software/miniforge3_x86
```

This now gives us two base conda environments.  One for installations native to the ARM architecture and native to M1 macs and one for all the other programs that haven't made an ARM version. The x86 programs will run using `Rosetta 2 translation environment`.

#### How to change between base conda installations
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

I put these in my home directory. To change between them all you have to do is source the file

```
source ~/condaX86
```

or

```
source ~/condaArm
```



## Qiime2

We do a lot of metagenomic analyses and Qiime can be one of the more challenging softwares to install.  In this tutorial, I will use it as an example of how to create a conda environment and install it on your local machine.  If you do not run metagenomic analyses yet then you can skip this setup.

* [Qiime install directions](https://docs.qiime2.org/2021.11/install/native/#miniconda)

```
source ~/condaX86

conda update conda
wget https://data.qiime2.org/distro/core/qiime2-2020.11-py36-osx-conda.yml
```

Here I am installing mamba for X86 as an install of mamba will initially run only for the ARM architecture, I had to create a specific environment to install the X86 version of mamba which I change into in order to use it to install X86 programs.  

The only downside to this is that I have to use the full path in order to source active the conda environment.
```
conda create --name mambaX86
conda install -n mambaX86 -c conda-forge mamba
mamba env create -n qiime2-2021.11 --file qiime2-2021.11-py38-osx-conda.yml
```

To run Qiime just activate the conda environment

```
# you have to use your full path in this case.
conda activate /Users/andrewseverin/GIFNew/Software/miniforge3_x86/envs/mambaX86/envs/qiime2-2021.11

```

[conda cheat sheet](https://kapeli.com/cheat_sheets/Conda.docset/Contents/Resources/Documents/index)

<!-- ## Other

* Singularity/Docker
* R packages
  * MetabolomicsR
  * tidyverse
* Conda environments?
* Bashrc
* eternal_history
* byobu
* qiime
  * which version am I currently using?
* camtasia
* juicebox
* latex Editor
* tassel
* Box Drive


```
.CFUserTextEncoding  .Rapp.history        .atom/               .zsh_sessions/       
.DS_Store            .Trash/              .zsh_history
``` -->


## Terminal setup
#### Giving Terminal Full Disk Access

If you start using Terminal, Mac will now ask you every time if you change directories permission to change into that directory.  To give Terminal access to all folders:

* Settings
  * Security & Privacy
    * Full Disk Access

#### Git files

You will want to copy over the following files from your old machine to your new machine to make github and Atom work again with pushes and pulls.  Github does not allow username/password pushes/fetches anymore.

* .gitconfig
* .git-credentials

I also had to create a new ssh key in my github account for my new laptop. If you need a refresher on how to do this see this [github tutorial](https://bioinformaticsworkbook.org/Appendix/github/introgithub.html#gsc.tab=0).


## .zshrc

```
# Aliases
pcat='pygmentize -f terminal256 -O style=native -g'
```
