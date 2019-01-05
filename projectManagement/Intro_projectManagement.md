# Introduction to Project Management

Project management is defined as the process of keeping track of the many aspects of a project from experimental design through publication. There are five main aspects to project management.

  * Folders and Files
  * Notebook record of commands, thought process and related code
  * Communication and Discussion
  * ToDo list
  * Record of your time spent

The overall goal of a solid project management process is the ability for a new person to join a project and be able to reproduce everything that has gone before and carry the project forward. This requires that all aspects of the project have the following properties.

* Reproducibility (version, hash, date of download)
* Location
* Date
* Version Control

With these properties in mind, let's step through each of the main aspects of project management.

## Folders and Files

Most Everyone uses folders to organize their files regardless of operating system (Windows, IOS, Linux).  A small change to how folders are labeled can significantly improve the ability for others to follow what you did, the order you did it and where you did it.

#### labeling

**Number folders in order of their creation by prepending the descriptive folder name that you usually create with a number starting with 00_, followed by 01_, 02_, 03_ and so forth like so:**

* Project folder
  * 00_RawData
  * 01_AnalysisStepOne
  * 02_AnalysisStepTwo
  * 03_AnalysisStepThree

The extra 0 from 00 - 09 keeps the folders ordered when executing the list (ls) command.  The you in 6 months or the next person on the project will thank you as they can numerically follow your steps based on the folder number and the date the folder was generated (ls -lha).

## Notebook

Science isn't reproducible if it isn't recorded.  Traditionally, paper notebooks have been used in the wet lab to record the science performed that day and the outcome.

In bioinformatics, science is recorded digitally in a digital notebook. We recommend using Github as a digital notebook.

Github is often thought of as a tool for programmers to version control code but it can also be used to version control text documents as an online notebook.  If you don't have a

If you are currently recording your notes by pasting commands in a text document or using R markdown or a ipython notebook then using a github private repo will be

We use [markdown](bioinformatics-workbook/Appendix/Markdown.md)









While there are many philosophies as to how to best organize your digital notebook, I will present one to you that works well.  

**For every folder you make, create a notebook that corresponds to that folder.**

* Project folder
  * 00_RawData.md
  * 01_AnalysisStepOne.md
  * 02_AnalysisStepTwo.md
  * 03_AnalysisStepThree.md



## Communication and Discussion



## ToDo Lists

## Record of time spent



Example use cases

* A grant was funded
* preliminary data generation for a grant proposal
*


## Learning Objectives



## Index

* Introduction to Slack
* Introduction to GitHub
*

# Introduction to using a github repo as an online notebook

The main points in reproducibility is for someone that has not worked on the project to know exactly what you did, when you did it, where you did it and the version you used.  This can be accomplished relatively easily by following these guidelines.

- Set up a folder on the computer or super computer that you will be analyzing your data.
- Create subfolders for each analysis step with a number associated with each folder starting at 00_.  This will allow someone to be able to follow the steps you took from the folder structure (examples below)
  - 00_rawData
  - 01_DataQC
  - 01_alignmentAndCounting
  - 02_RNA-SeqDE
- Create a corresponding Mark down file in this github repository that corresponds to the work found in each folder.
- For every step you make record it in this mark down folder along with
  - title
  - date
  - Computername and path of working directory
