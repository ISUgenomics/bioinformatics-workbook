# Guide for installing various types of programs in Linux

This handy guide is for installing programs in `UNIX` environment. Most of these steps assume that you are installing package in a group accessible location, without root access and utilizing the environment module systems for package management. However, you can easily modify these steps for other cases as well.



### Conda

One of the easiest ways you can install you own software in your home or project directory is through the [Conda package manager](https://conda.io/docs/user-guide/getting-started.html). Thousands of biological packages and their dependencies can be installed with a single command using the [Bioconda repository](https://bioconda.github.io/) for the Conda package manager.


### Unpacking

Packages are usually compressed in many different ways for easy handling. Before proceeding to installation, it must be unpacked. Depending on the compression (looking at the extension) use any of the following commands to unpack.
Although `tar` can auto detect the compression type and decompress the archive with the `-xf` options, you can also specify what type compressed files you're providing. For most cases `tar -xf` will do the trick:

```
tar xf package.tar.gz
```

But to be specific you can:

```
# Files having *.tar.bz2 extension
tar xvjf package.tar.bz2
# Files having *.tar.gz extension
tar xvzf package.tar.gz
# Files having *.bz2 extension
bunzip2 package.bz2
# Files having *.rar extension
unrar x package.rar
# Files having *.gz extension
gunzip package.gz
# Files having *.tar extension
tar xvf package.tar
# Files having *.tbz2 extension
tar xvjf package.tbz2
# Files having *.tgz extension
tar xvzf package.tgz
# Files having *.zip extension
unzip package.zip
# Files having *.Z extension
uncompress package.Z
# Files having *.7z extension
7z x package.7z
```
Once de-compressed, proceed with installation, depending on what type of package you are installing.

### Regular Linux packages

#### Standard approach

Many `Linux/UNIX` programs comes with a standard set of files that lets you install programs with ease. After unpacking, if you see `configure` file in unpacked directory, use this approach.
```
./configure --prefix=/group/accessible/location/packagename
# once complete you'll see 'Makefile' in that directory
make # reads instructions from 'Makefile' and builds executable programs
make check # not needed, but helps in troubleshooting
make install # installs program using 'Makefile' directions
```
In case if something goes wrong or you get an error saying that you need `package x` before installing, then you can undo these steps before attempting installation again.
```
make clean
```
If the program doesn't work as intended or something goes wrong after installation, many programs can be safely uninstalled
```
make uninstall
```
Note that all of these options are not supported by all makefiles so your mileage may vary. It is also good idea to run all the above command in a `build` directroy inside the package directroy, so that if something doesn't work you can easily delete the `build` directory to start over.


#### Special cases

*1. No configure file*

Some programs will already have a `Makefile`. These programs do not need the first step (running `configure`), you can simply install them by typing
```
make
make install
```
The executables are generally created either in the same directory or in the `bin` directory, within the package directory. Sometime these packages will allow you to install other locations as well. Consult the `README` or `INSTALL` files that came with the program or edit the `Makefile` to hard code the installation directory. In some cases, setting `PREFIX` variable to the desired installation location will also do the trick.

```
PREFIX=/your/installation/location  make
```


*2. Source obtained as RPM package*

In some rare cases when you don't find a package for Red Hat Linux, but you have `rpm` package for the program, then use these steps to extract and install:
  * Find the correct RPM package for your system. This [[http://pkgs.org | webiste]] lists all RPMs available, and are free to download. All CentOS RPM's work on Red Hat. Download the `Source` package (not Binary).
  * Extract the package: use `rpm2cpio package.rpm |cpio -idv` command to extract, you should see a `*tar.gz` or other type of compressed program, if this completes successfully.
  * Follow the steps described earlier to install like any other package.
  * In rare cases, when you have patches (extracted from RPM), you might have to apply it before you install. After extracting `*.tar.gz` file, `cd` in to the directory and run `patch -Np1 -i path/to/file.patch` and install as usual.

*3. Needs cmake*

If the README file says that you need to use `camke` command, then use these steps to install:
```
# after extraction, cd to the package
cd package
mkdir build
cd build
cmake ..
# if you want it in a different directory, then
cmake -DCMAKE_INSTALL_PREFIX:PATH=/location/for/installation ..
make
# if this completes successfully, you will see a bin folder above this current directory
# that will have the executables
```

## Group installation
Rather than trying to install to root, it is useful to be able to install programs that are common for your group/lab.  To attain this type of installation it is necessary to download and load your own version of `python` (not the `python` installed in root, by admin). Global installations are not possible when you are not the administrator, and personal directory installations are only available to one person.  When it is desirable to install to group.  Start with Python installation in your groups directory (example: /work/GIF on condo).

### Python packages
Using our own `python` will allow writing/installing modules to it as needed. After unpacking, `cd` to the package, and install it as follows:
```
module load python
python setup.py install # all executables will be stored in python/bin (not in package directory)
```
If in case if you need to test out something and not install it as module, you can install in a personal location as well:
```
python setup.py install --local=/home/username/mydir
# or simply as
python setup.py install --user # executable's will be in ~/.local/bin directory
```

Any package available at [[https://pypi.python.org/pypi | PyPi ]] can be managed using these commands as well
```
module load python
pip install SomePackage # installs a python package
pip show --files SomePackage # shows what files are installed for the particular package
pip list --outdated # lists what packages are outdated
pip install --upgrade SomePackage # upgrades a package
pip uninstall SomePackage # uninstalls a package
pip freeze # lists all the packages that are currently installed and their version
```

### Perl modules
Once the module is loaded, use the following set of commands to install any `perl` modules.
```
module load perl
```
If there is a `Makefile.PL`
```
 perl Makefile.PL PREFIX=/home/users/dag   # makes the system specific makefile
 make          # builds all the libaries
 make test     # runs a short test
 make install  # installs the package correctly.
```
If there is a `Build.PL`
```
perl Buil.PL
./Build test
./Build install
```
The module will be installed in the group's perl folder (not in the package directory). So, like you did in `Python` you need to set up a dummy module file that load `Perl`

### R or Bioconductor packages

Installing `R` libraries for the group is really easy since you don't have to do anything different from the way you install packages to your home directory. GIF has its own `R` version installed as module and it is configured such that it will automatically install the package in the correct location, when you are using this module.
```
module load R
R
# R command prompt will appear
>
```

#### Installing CRAN R Packages
CRAN packages are by far the easiest. From within R prompt, type:
```
module load R
R
# R command prompt will appear
> install.packages("some_package") # include quotes!
# if it prompts to select the closest mirror, choose IA, which is `77`
```
Once installed, you will be back at `R` prompt, load the installed package to see is everything is fine.
```
library(some_package)
# this should load the package and return without any error message
```
#### Install manually downloaded R Package
Some packages that aren't in CRAN but are available from the author directly, can be installed for group as well. Download the the `package.tar.gz` from the author's website.
```
module load R
R CMD INSTALL package.tar.gz
# This will install the package for the group.
# Check to see if it works
R
# R command prompt will appear
> library(package)
# this should load the package and return without any error message
```

#### Installing Bioconductor Modules
For Bioconductor packages, follow these steps:
```
module load R
R
# R command prompt will appear
> source("http://www.bioconductor.org/biocLite.R")
> biocLite(c("package_name"), dependencies=TRUE)
# check if package is installed
> library(package_name)
# this should load the package and return without any error message
```

#### Managing packages

*List available packages*

To get a complete list of packages that are already installed, load the `R` module and enter the R prompt. From there, type the following command:
```
> library()
```
To get all packages installed along with their version number, type
```
> installed.packages()[,c("Package","Version")]
```

#### Upgrade packages
For the R packages that you installed from `CRAN` can all be upgrades in single command
```
update.packages() # upgrades all packages
package.status() # says if 'ok' (no updates), 'upgrade' (needs update) or 'unavailable' (package removed from repository)
```
Other useful option to check the status of all packages currently installed is:
```
> inst <- packageStatus()$inst
> inst[inst$Status != "ok", c("Package", "Version", "Status")]
```
#### Uninstall package

Packages can be uninstalled easily using `remove.packages` command
```
> remove.packages("package_name")
```

### Java programs
Precompiled java programs that come as `.jar` files, can be placed in any directory and can be called from there. For using it with environment modulefile, you need to do these steps. First, create directory (program name) and sub-directory (version number). Place the `.jar` file in this sub-directory. Within this create another directory and call it as bin. For all `.jar` files in `/programname/version/` create a text file in `/programname/version/bin`. This text file will just have a single line, something like:
```
java program_name.jar
```
Change permission for these text files so that they can be executed.
```
chmod +x -R /programname/version/bin
```
In your module file, you need to add this line:
```
prepend-path PATH /programname/version/bin
```
Now the `.jar` files can be simply called as `'programname` (once module is loaded). No need to add `java` in front.
