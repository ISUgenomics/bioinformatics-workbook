---
title: "Project Management"
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Project Management


# Introduction to Project Management

Project management is an ongoing learning endeavor to discover what works well for you and your team.  What is presented here is a project management system that has evolved in the Genome Informatics Facility at Iowa State University, which is tasked with analyzing dozens of projects every year.  

Project management is defined as the process of keeping track of the many inputs and outputs of a project from experimental design through publication.

Many components of a project are generated, in an attempt to organize the inputs and outputs of a project.  How these components relate to each other and their organized is directly related to project management.

The main components of project management include but may not be limited to:

  1. Folders and Files
  2. Notebook record of commands, thought process and related code
  3. Backup of raw data and analysis folders
  4. Communication and Discussion
  5. ToDo list
  6. Record of your time spent


The rest of this tutorial steps through the components of project management, how they are related and importantly how the adoption of some very simple rules help integrate these components and do not come with a large learning curve that can be prohibitive to quick adoption by your team.

## 1. Folders and Files

Most people use folders to organize their files regardless of operating system (Windows, IOS, Linux).  A small change to how folders are labeled can significantly improve the ability for others to follow what you did, the order you did it in and where you did it.


**Number folders in order of their creation**

Prepend the descriptive folder name that you usually create with a number starting with 00_, followed by 01_, 02_, 03_ and so forth like so:

* Project folder
  * 00_RawData
  * 01_AnalysisStepOne
  * 02_AnalysisStepTwo
  * 03_AnalysisStepThree

The extra 0 from 00 - 09 keeps the folders in numerical order when executing the list ```ls``` command up to 99 folders.  The you in 6 months or the next person on the project will thank you as it is easy to numerically follow the steps based on the folder number and the date the folder was generated ```ls -lha```.

The name of your project folder should be descriptive of the project, for example, 2020_Severin_UnicornRNA-Seq or 2020_Severin_SeriolaGWAS.  Since I work with many different researchers, I tend to name my Github Repos and project folders with the Year_ResearcherName_shortProjectDescription.

## 2. Notebook

Science isn't reproducible if it isn't recorded.  Traditionally, paper notebooks have been used in the wet lab to record the science performed that day and the outcome.
In bioinformatics, science is recorded digitally in a digital notebook. We recommend using GitHub as your digital notebook.

#### GitHub

GitHub is often thought of as a tool for programmers to version control code but it can also be used to version control text documents as an online notebook.  If you don't have a GitHub account, you can set one up based on this [GitHub Tutorial](../Appendix/github/introgithub.md) and the free version of GitHub now allows for unlimited free repositories. In addition, as a researcher or educator you can request a free upgrade to get access to [unlimited users in private repos.](https://help.github.com/articles/about-github-education-for-educators-and-researchers/).  

One of the great aspects of GitHub is that you can use [markdown](../Appendix/Markdown.md) syntax to very quickly format your notebook.  For example, This tutorial is written in markdown and you can make headers```# ## ###```, bullet points ```*```, code blocks and insert images very easily.  Markdown can also be easily converted to PDF for sharing.  It is also straightforward enough to share the repo itself when passing on a project or collaborating with a student or colleague.

Assuming you now have a GitHub account, create a GitHub repository with the same name you labeled your project folder (ie: 2020_Severin_UnicornRNA-Seq or 2020_Severin_SeriolaGWAS).  

<span style="color:Green">Here is the first integration between project components **(Folder names match file names in your GitHub notebook)**.  This is where the simplicity of this approach gives dividends and returns on such a small investment in time.</span>

**For every subfolder create a markdown file**

For every subfolder you make in your project folder on your remote machine or computer, create a markdown file in your GitHub repository.

* **Github Repo**
  * Notebook_LastName
    * 00_RawData.md
    * 01_AnalysisStepOne.md
    * 02_AnalysisStepTwo.md
    * 03_AnalysisStepThree.md

After you git clone your repo you can create a folder named Notebook_LastName and then add markdown files corresponding to each new subfolder you add during the analysis pipeline.

If you work on multiple machines and your laptop then create a folder for each machine.

* Notebook_LastName
  * remote machine
    * 00_RawData.md
    * 01_AnalysisStepOne.md
    * 02_AnalysisStepTwo.md
    * 03_AnalysisStepThree.md
  * laptop
    * 01_AnalysisStepOne.md
    * 02_AnalysisStepTwo.md

<span style="color:Blue">Why not just put the markdown files on the main level?</span>
  * It is cleaner and there are other files we want at the main level
  * More than one person may be working on this project and we want all those files and folders to be separated into the User's notebook.

#### Example Repo
For every new project, I copy over this starter repo [Repo_skeleton](https://github.com/ISUgenomics/Repo_skeleton).  It has some other markdown file suggestions that include.


| File | Description|
| :-- | :-- |
| 00_Files.md | a list of all important files and where they are located|
|00a_MetaData.md | A description of the metadata for each sample|
|01_Background.md | a place to put more information about the organism and the project objectives|
| 02_Methods.md | Write up your methods as you go|
| 03_Results.md | Place the results here as you go|
| 04_Introduction.md | Once the Methods and Results are written up start collaborating on an Introduction|
| 05_Discussion.md | Write a rough draft of the Discussion or put ideas here for later writing of the discussion|
| 06_AuthorInfo.md | Good to be thinking about authorship early|
| Notes.md | Other Miscellaneous Notes|


#### Atom (text editor)  <- Download this tool

Using the GitHub interface on the web or even a basic text editor is not very efficient. I strongly recommend downloading and installing the [GitHub ATOM](https://atom.io/) markdown editor.  This is a text editor that understands markdown syntax and allows you to preview markdown that you write in an adjacent window.  See [markdown](../Appendix/Markdown.md) tutorial for more information.


#### What, When, Where, Why, How and the Result
For every step you perform bioinformatically, record it in your markdown file along with
  - Title
  - Date
  - Computername and path of working directory
  - Description of what you are trying to do and why
  - Code block of the commands you used
  - Description of the result

This last step is important as it is a discrete block of notes that is the smallest unit of replicable notes in your notebook.  The Why, How and Result can be replicated multiple times if it is the same day in the same folder.

#### <span style="color:Green">Example</span>
Inside my notebook folder inside my computer folder in a markdown file called 01_QC.md that corresponds to my 01_QC folder where I am performing a Quality Control of raw data I might find. the following

* notebook_Severin
  * Condo
    * 01_QC.md

````
# Quality check on Unicorn RNASeq data             (What)

* Jan 19, 2021                                     (When)
* /work/GIF/severin/Unicorn_RNASeq/01_QC           (Where)
                                                   (Why)
I want to make sure that the quality of my reads obtained from the Illumina Novaseq look good before proceeding.   


## FASTQC                                          (How)

```
module load fastqc/0.11.7-d5mgqc7
module load parallel/20170322-36gxsog

mkdir fastqcOutput
  parallel "fastqc {} -o fastqcOutput" ::: ../00_rawdata/*.fastq.gz
```

## MultiQC                                         (How)

```
#load python 3.6
module load python/3.6.3-u4oaxsb
# multiqc wasn't working so I reinstalled it
pip uninstall multiqc
pip install multiqc
multiqc --version
multiqc, version 1.8
multiqc .
```
                                                    (Result)
The output from FASTQC and multiqc look good!  Proceeding to differential expression analysis of unicorn horn between activated and unactivated samples.

````
As you can see from the example above, when we implement the what, when, where, why, how and result note block, it increases comprehension for anyone including yourself 6 months from now or the new researcher taking over your project.  This note block should be repeated for every analysis, analysis attempt, Error trouble shooting, parameter changes.  For any (How) section that advanced the analysis or may be useful in the future, surround that how with the (What, When, Where,Why)(How)(Result/Discussion).

**Note:** <span style="color:Blue">Always record the version of the software you are using.  All publications now require this information for reproducibility and it can be a pain to find after 6 months and several potential version updates of said software. </span>

I will admit that I am not perfect at this and sometimes forget to add the when or where but I have discovered the closer I follow this note block rubric the easier it is to find information later.

## 3. Analysis folder backup and cleanup

Backup of your raw data and analyses folders is an important part of the life cycle of a project.  When you have as many projects as I see over the course of the year, the amount of storage required tends towards 100 Terabytes.  If my group does not stay on top of backing up and returning completed projects, the available space to perform new projects becomes a limiting factor. Below you will find our general process.

1. After completing a project and sending the data back to the researcher. The folder is cleaned up of any analyses that are safe to delete.
  a. Failed attempt folders that are no longer needed.
  b. Alignment intermediates (sam) before bam file generation
  c. Split files that may have been generated before combining to generate a final combined output.
  d. and similar types of files or outputs that can get generated.  Watchout for programs like Trinity as they can generate hundreds of thousands<span style="color:Red"> (300,000-1,000,000) </span>files.  Those files should be removed (cleaned up).
2. Raw data should be backed up to a secondary site.
  In our case, we have access to cloud services such as Box and Google drive that allow us to backup data and folders.  We also have a Large Scale Storage solution through our University that allows us to purchase affordable hardrive space.  Where you back up your projects will be dependent on your location.  However, it is good practice <span style="color:Red">not</span> to store your raw data and final results on the HPC resource as most HPC resources do not have redundant backup and are replaced every 3-5 years.

3. After 3-6 months or publication of the paper, whichever comes first, the project folders that you haven't touched should be further cleaned with additional intermediates removed. These may include:
    * alignment files that are no longer needed because you have the vcf file
    * other files depending on the type of analysis.

It is good to keep the general structure of the analysis and files such as logs and submission scripts.

This process a little thought if you are using cloud storage because many cloud storage solutions have limits on the number of files you can upload to a shared folder.  Google for example has a limit of 400,000 files and the file size has to be less than 5 Terabytes. There is also a daily upload limit of 768 Gb/day which will continue to try until it is uploaded by will slow down your transfer limiting you to the 768 Gb/day.

#### Useful Scripts for folder cleanup

Here are a few helpful commands for getting a folder in shape for saving to a secondary location.

* filecount.sh

This bash script will generate a count of the number of files in each subfolder.
```
#!/bin/bash

find . -maxdepth 1 -type d | while read -r dir; do printf "%s:\t" "$dir"; find "$dir" -type f | wc -l; done > fc.txt &  
```

Any folders that have more than 1000 files should be gzipped before rcloning to a secondary location and finally folders should be removed from the HPC system

Example output

```
more fc.txt
.:      45795
./ProjectA:       93
./scratch:      62
./OrthoReD:     122
./scripts:      6
./ProjectB:      368
./Unicorn_RNASeq: 24003
./L3:   7619
./mutantDNA:       166
./SuperPowerStudy: 4226
./student:    5722
./facts:     82
./nr:       3322
```

* File Size

File sizes for all files and folders in your project folder can be determined with this command

`du -hs * > fileSizes.txt`

* File list

It is often helpful to create a text file that contains all the files that are in the folder.  This is especially true if you end up zipping many of the subfolders.  If you want to look for a file later, you can first search the file list and see if it is there rather than having to download the zipped folder, unzipping and then searching for the file you want.

`find . > files.txt`

#### Summary

These three commands are a great help for cleaning up directories

```
filecount.sh
du -hs * > fileSizes.txt
find . > files.txt
```

**NOTE:** filecount.sh will need to be in your PATH in order for you to just call it as I do above.

## 4. Communication and Discussion

In a group with multiple active projects, it is helpful to have a running dialogue of progress and discussion for each project that can be easily referred to and searched.  Many times this comes in the form of email exchanges.  

#### Slack

However, a communication platform based around group chats can also be extremely effective when you have a medium sized group (3+) of knowledgeable individuals.  We use

[Slack](../Appendix/slack.md) and create a channel for each active project.  This group chat (text messaging style) communication generates a transcript of conversations about projects that you can go back and reread. It also permits, text, images and files to be attached during conversation.  I strongly encourage you to use it primarily for text communication and discussion and not for file storage.

You can create your own team using [this link: Team signup](https://slack.com/r/0mklqxar-0tnh0jhj) and get a $100 credit if you decide to purchase it beyond the free version.  The upgraded version will let you search beyond the 10,000 most recent messages in your group and there is a steep academic discount.

#### OneNote

Microsoft Onenote has the ability to create notebooks, tabs and pages. This can be done on your laptop or tablet using a pen.  Both will result in a significant reduction in your paper usage.  In addition you will have a well organized record of every in person meeting you have with a student or colleague on a particular project.

* Meetings Notebook
  * ProjectName Tab
      * Meeting date1 page
      * Meeting date2 page
      * Meeting date3 page
  * ProjectName2 Tab
    * Meeting date1 page
    * Meeting date2 page
    * Meeting date3 page

Since the Genome Informatics Facility works on multiple projects with multiple PIs, we use the PI name for the project name tab.  Within a lab group you may have a folder for every person in your group.

## 4. ToDo Lists

I used to write todo lists on a piece of paper and continue to add more to the todo list each day and cross off the items I completed.  After a while, this todo list was over a page long and important items would get buried far down on the list.  In order to change priority of the items on the list I had to rewrite the entire list or create a new sublist on a separate piece of paper.  This was frustrating and time consuming.  There is a better way.

Use a Kanban board. Kanban boards can be thought of as a 2D ToDo list. Take a whiteboard and section it into 6 columns. (There is a digital equivalent that I will provide you below so don't actually do this with a real whiteboard)

* New tasks
* IceBox
* Backlog
* In progress
* In review
* Completed

Take a post-it note and write down the task you need to do and everything you need to know about how to do it.  Repeat for every task you have. Place all post-it notes into the ```New tasks``` column on your whiteboard.

Now that all the tasks are on the whiteboard, separate out those tasks that are higher priority from those that are lower priority and place them into ```In Progress``` and ```Backlog```, respectively.

The ```Icebox``` column is for those ideas that are rumbling around in your head that are really cool and you don't want to forget but don't have time for right now.  Write down that idea/task and put it in the ```Icebox``` column for when you do have time and get it out of your head.

Now prioritize the tasks in ```In Progress``` and ```Backlog``` columns by placing the most important tasks at the top of the column.

When you complete a task move it to ```In review``` or ```Completed``` columns.

Cool idea right?  But a real-life whiteboard and post-it notes would be very time consuming so we use a digital whiteboard.

[Zenhub](https://www.zenhub.com/pricing) is free for public, personal, and academic repositories and the chrome or firefox plugin can be downloaded [from here](https://www.zenhub.com/extension).  It uses the issues in a github repository.  So now you can for every project have a Kanban board for all tasks related to that project.  If you have more than one project it may be helpful to create a separate private repository just for your todo lists on several projects.


## 5. Record of time spent

Bioinformatics projects can be time consuming and frustrating with a lot of trial and error.  It often happens that at the end of a week, it doesn't feel like anything was accomplished.  Keeping track of your time spent on a daily basis can bring peace of mind.  This can be easily accomplished using a Google spreadsheet with days of the month at the top as columns and project names as rows on the left. For each project there should be two rows; hours spent and description of what you did.  The column widths can be relatively short since hours spent in a day is not more than 2 digits and the description can just flow into adjacent cells.  It is also helpful to see where your time is being spent and better manage your work/life balance and provide a tool for those in your group to manage their work/life balance.

Here is an example of a [google sheet](https://docs.google.com/spreadsheets/d/1Lj2xXjnWkx8ww7l23OjRET77J3_jq6Cj61Iu6bebMxw/edit#gid=0) used for keep track of time.
