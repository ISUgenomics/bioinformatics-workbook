# Macbook Pro Installation

## Optional Software
### Word Processing
* Microsoft Office Suite (Word, Excel, Powerpoint, Teams, OneNote)
* Apple Suite (Keynote, Numbers, Pages)

### Password Manager
* Password Manager (Dashlane, OnePass, etc)

### Screen Capture/Video Editors
* [Camtasia](https://www.techsmith.com/download/camtasia/)
* [Apple Pro Apps bundle for education](https://www.apple.com/us-edu/shop/product/BMGE2Z/A/pro-apps-bundle-for-education)
  * Final Cut pro
  * Compressor
  * Motion
  * Logic Pro

## Team communication

* [Slack](https://slack.com/downloads/mac)
  * The easiest way to reload your teams is to sign in on the slack website and then use your email to sign in for each team you want to be signed in on your desktop application.
* [Zoom](https://zoom.us/support/download?os=mac)
  * Don't forget to copy over your zoom backgrounds



**
## Install Xcode Developer tools

In order to get the most out of your macbook, this developer's toolkit is required.  Installing XCode at a minimum also installs git for github commands.

* The download takes a long time
* The unzipping takes a long time
* Move the Xcode App to your applications

Install the XCode Command line tools

```
xcode-select --install
```

## Install Rosetta 2
For those applications that are still compiled for Intel processors, Apple has Rosetta 2

Run the following command on the terminal

```
/usr/sbin/softwareupdate --install-rosetta --agree-to-license
```

Restart your computer

## Install Atom editor

* https://atom.io

Atom did not open for me the first time when I double clicked on it.

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

If you need a vpn to access your HPC you will need to download that.  For us at Iowa State we use the Cisco Any Connect client

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

* OneTab (you can download this from App extensions store from Apple)

## Install R

Install R from the package from R cran site.

* [Install R](https://mirror.las.iastate.edu/CRAN/)
  * https://mirror.las.iastate.edu/CRAN/
    * https://mirror.las.iastate.edu/CRAN/bin/macosx/big-sur-arm64/base/R-4.1.2-arm64.pkg

I got the following error when I tried to install from DMG so don't do that.
**Error** "R can't be opened because Apple cannot check it for malicious software"
  * [Monterey Arm version 4.2](https://mac.r-project.org/monterey/R-devel/R-GUI-8008-4.2-monterey-arm64-Release.dmg)

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

* https://cran.r-project.org/web/packages/renv/
* https://6chaoran.wordpress.com/2020/07/20/introduction-of-renv-package/


pacman?
R environment manager

mamba conda

update

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

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> /Users/YourNAME/.zprofile
   eval "$(/opt/homebrew/bin/brew shellenv)"

```

## Conda

Install Miniforge3 that emphasizes the arm architecture.

* [Miniforge3](https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh)

## Qiime2

* [Qiime install directions](https://docs.qiime2.org/2021.11/install/native/#miniconda)

```
conda update conda
```



## Other

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
* figtree
* igv
* juicebox
* jbrowse desktop
* latex Editor
* pymol
* tassel

* Box Drive


```
.CFUserTextEncoding  .Rapp.history        .atom/               .zsh_sessions/       
.DS_Store            .Trash/              .zsh_history
```


## Terminal setup
#### Giving Terminal Full Disk Access

If you start using Terminal, Mac will now ask you every time if you change directories permission to change into that directory.  To give Terminal access to all folders:

* Settings
  * Security & Privacy
    * Full Disk Access
