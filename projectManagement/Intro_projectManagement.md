---
title: "Project Management"
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

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


**Number folders in order of their creation by prepending the descriptive folder name that you usually create with a number starting with 00_, followed by 01_, 02_, 03_ and so forth like so:**

* Project folder
  * 00_RawData
  * 01_AnalysisStepOne
  * 02_AnalysisStepTwo
  * 03_AnalysisStepThree

The extra 0 from 00 - 09 keeps the folders ordered when executing the list (ls) command.  The you in 6 months or the next person on the project will thank you as it is easy to numerically follow the steps based on the folder number and the date the folder was generated (ls -lha).

The name of your project folder should be descriptive of the project, for example, WhiteAbaloneRNA-Seq or SeriolaGWAS.

## Notebook

Science isn't reproducible if it isn't recorded.  Traditionally, paper notebooks have been used in the wet lab to record the science performed that day and the outcome.

In bioinformatics, science is recorded digitally in a digital notebook. We recommend using Github as a digital notebook.

#### Github

Github is often thought of as a tool for programmers to version control code but it can also be used to version control text documents as an online notebook.  If you don't have a Github account, you can set one up based on this [Github Tutorial](/Appendix/github/introgithub.md) and as a researcher or educator you can request a free upgrade to get access to [unlimited private repos.](https://help.github.com/articles/about-github-education-for-educators-and-researchers/)

One of the great aspects of Github is that you can use [markdown](/Appendix/Markdown.md) syntax to very quickly format your notebook.  For example, This tutorial is written in markdown and you can make headers, bullet points, code blocks and insert images very easily.

Assuming you now have a github account, create a github repository with the same name you labeled your project folder (ie: WhiteAbaloneRNA-Seq or SeriolaGWAS).  

**For every subfolder you make in your project folder on your remote machine, create a markdown file in your github repository.**

* Notebook_LastName
    * 00_RawData.md
    * 01_AnalysisStepOne.md
    * 02_AnalysisStepTwo.md
    * 03_AnalysisStepThree.md

After you git clone your repo you can create a folder named Notebook_LastName and then add markdown files corresponding to each new subfolder you add during the analysis pipeline.

#### Atom

If you haven't already, download and install [Github ATOM](https://atom.io/).  This is a text editor that understands markdown syntax that allows you to preview markdown that you write in an adjacent window.  See [markdown](/Appendix/Markdown.md) tutorial for more information.

#### What, When, Where, Why, How and the Result
For every step you make record it in your markdown folder along with
  - Title
  - Date
  - Computername and path of working directory
  - Description of what you are trying to do and why
  - Code block of the commands you used
  - Description of the result

## Communication and Discussion

In a group with multiple active projects, it is helpful to have a running dialogue of progress and discussion for each project that can be easily referred to and searched.  

Slack is a communication platform that is based around group chat channels. It is a transcript of conversations about projects that you can go back and reread. Slack permits, text, images and files to be attached during conversation.


## ToDo Lists



## Record of time spent



Example use cases

* A grant was funded
* preliminary data generation for a grant proposal



## Learning Objectives



## Index

* Introduction to Slack
* Introduction to GitHub
* Introduction to Markdown
