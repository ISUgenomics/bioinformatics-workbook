---
title: "Useful Programs and Unix Basics"
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

## Creating your own Singularity containers.  What you need to know.


This tutorial will help you create your own Singularity container using a github, a recipe file and Singularity Hub.  It will also suggest a GitHub repo organization that will maximize a researcher's ability to utilize your scripts and pipelines in your repo.

You will need:
* A GitHub Repo [https://github.com/](https://github.com/)
* A Singularity Account [https://www.singularity-hub.org/](https://www.singularity-hub.org/)
* Singularity installed locally via the vagrant virtual machine or on a linux machine where you have root access. [See Intro to Singularity](Intro_Singularity.md)


#### Example recipe (shub://ISUGIFsingularity/utilities:utilities.1.0.0) for some basic scripts we use in GIF
While it isn't necessary to build a container for this as it is more efficient to have the scripts in your path on your machine, I used this a toy example to familiarize myself with container images, recipes and wrappers.

Singularity Recipes have Five main sections
* **Header information**
* **%labels**       :information about the container
* **%environment**  : environmental variables and modules to load.
* **%post**         : commands to run during the creation of the container
* **%runscript**    : program to run if the image is executed        

You can find more information about these sections on the [singularity website](http://singularity.lbl.gov/docs-recipes#environment)

You can bootstrap from existing images either from Docker or Singularity.
```
Bootstrap:docker
From:centos:7

Bootstrap:shub
From:ResearchIT/spack-singularity:openmpi
```

Here is an example where we use a prebuilt singularity container that has spack installed.  [Spack](https://spack.io/) is community driven package manager that attempts to improve scientific reproducibility by providing installation methods for software installation. By combining spack methods with container recipes we can capture both the reproducibility of installation but also the entire environment.  In addition, spack makes creating singularity containers extremely straight forward so long as the [package is found in spack](http://spack.readthedocs.io/en/latest/package_list.html).

#### Example prepare_genome_modules

In this example we use the spack base singularity image and build onto it by installing gsnap, bowtie2, bwa, gatk bedtools2, samtools, picard and java.  These installations which can be long and involved are a 3 word command in spack (see %post section below)

```
Bootstrap:shub
From:ResearchIT/spack-singularity:openmpi

%labels
MAINTAINER severin@iastate.edu
APPLICATION genomeModules

%help
This container contains all the necessary programs to create genome modules.
See https://github.com/ISUGIFsingularity/genomeModules.git for more information

%environment
source /etc/profile.d/modules.sh
SPACK_ROOT=/opt/spack
export SPACK_ROOT
export PATH=$SPACK_ROOT/bin:$PATH
source /etc/profile.d/modules.sh
source $SPACK_ROOT/share/spack/setup-env.sh
export PATH=$SPACK_ROOT/isugif/utilities/bin:$SPACK_ROOT/utilities/wrappers:$PATH
#for d in /opt/spack/opt/spack/linux-centos7-x86_64/gcc-4.8.5/*/bin; do export PATH="$PATH:$d"; done

module load gmap-gsnap
module load bowtie2
module load bwa
module load gatk
module load bedtools2
module load samtools
module load picard
module load jdk
export _JAVA_OPTIONS="-Xmx100G"
module load parallel

%post
export SPACK_ROOT=/opt/spack
export SPACK_ROOT
export PATH=$SPACK_ROOT/bin:$PATH

yum -y install bc paste
yum clean all

export FORCE_UNSAFE_CONFIGURE=1

source $SPACK_ROOT/share/spack/setup-env.sh
spack install picard
spack install gmap-gsnap
spack install bowtie2
spack install bwa
spack install gatk
spack install bedtools2
spack install samtools
spack install parallel

#for d in /opt/spack/opt/spack/linux-centos7-x86_64/gcc-4.8.5/*/bin; do export PATH="$PATH:$d"; done


cd $SPACK_ROOT

%runscript
echo "This container contains a environment and all prequisite programs to run prepare_genome_modules.sh"


```


#### Building an image
If working locally start a singularity VM

```
vagrant init singularityware/singularity-2.4
vagrant up
vagrant ssh
```

Place the example recipe above into a file named recipe then create a singularity image we will call test.simg using the following command.

```
sudo singularity build test.simg recipe
```


## Singularity containers and GitHub Repos

One of the goals of GitHub is to share code and scripts that we find useful so that others may also find it useful.  Containers take this to another level by providing all the prerequisites built into the script itself.  ISUGIF is implementing the use of containers for all our scripts we have online so that they can be used without fear of problems due to differences in HPC environment like different names for common required modules.  You can find our [singularity repos here](https://github.com/ISUGIFsingularity)

#### Organization of the GitHub Repo

let's start by downloading the example we have used throughout this tutorial.

```
git clone https://github.com/ISUGIFsingularity/utilities.git
cd utilities
ls
```
```
README.md         Singularity.1.0.0 utilities         wrappers
```

In this repo, you will find the
* **README.md** Text file that describes the repo which also includes where to download the singularity from singulairty-hub
* **Singularity.1.0.1** Text file with the recipe that Singularity hub uses to generate the singularity container on their website
* **utilities** A folder that contains all the scripts I wanted to share
* **wrappers** A folder that contains bash script wrappers that call the singularity program using the singularity image produced by Singularity-hub and the scripts in the utilities folder.

The GitHub repo can now be used in the traditional way where you git clone the repo and use the functions inside utilities or it can be used with the new method by calling singularity containers.  The new method has the advantage that all modules and prerequisites you may need in a script or pipeline is already included!  

#### creating wrappers for the singularity commands

Here is an example bash script wrapper for a singularity execution of a command in a container.  The use of singularity is completely hidden to the user.

```
new_assemblathon
#!/bin/bash

singularity exec ISUGIFsingularity-utilities-master.simg ../utilities/new_Assemblathon.pl
```


#### Naming Schema

* Singularity requires that you have a recepe named as Singularity.XXX  Where XXX is usually the version number.  
* These Recipe files are contained in individual github Repos with more descriptive names.
* On Singularityhub it will be OrganizationName/RepoName/recipeFilename
* On Github you will create a descriptive named repo for every singularity container you want to make and create a recipe file contained in that Repo starting with Singularity.


#### More Information about containers

* [Apps within containers](http://singularity.lbl.gov/docs-scif-apps)
* [Recipe sections](http://singularity.lbl.gov/docs-recipes#environment)


#### Next tutorial Modifying existing containers for interesting

* [Modifying Existing Containers](modifyingExistingContainers.md)

---
[Table of contents](../../programs.md)
