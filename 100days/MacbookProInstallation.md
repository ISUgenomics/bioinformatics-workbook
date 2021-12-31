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
 



** 
## Install Xcode Developer tools

In order to get the most out of your macbook, this developer's toolkit is required.  Installing XCode at a minimum also installs git for github commands. 

* The download takes a long time
* The unzipping takes a long time
* Move the Xcode App to your applications

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


## Install R

Install R from the package from R cran site.

* [Install R](https://mirror.las.iastate.edu/CRAN/)
  * https://mirror.las.iastate.edu/CRAN/
    * https://mirror.las.iastate.edu/CRAN/bin/macosx/big-sur-arm64/base/R-4.1.2-arm64.pkg

I got the following error when I tried to install from DMG so don't do that.
**Error** "R can't be opened because Apple cannot check it for malicious software"
  * [Monterey Arm version 4.2](https://mac.r-project.org/monterey/R-devel/R-GUI-8008-4.2-monterey-arm64-Release.dmg)

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
