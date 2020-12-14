---
title: "Workflow tutorials index"
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Workflows

Workflow languages are becoming common practice in bioinformatics and many big data sciences to automate pipelines.  As a person advances in your bioinformatics knowledge, they tend to create many scripts commonly in bash, python, perl, insert your favorite language, etc to make the processing of data easier.  Longer scripts can be created in these languages but each project is unique and tend to have exceptions that would have to be accounted for resulting in more and more complicated script.  Workflows permit the user to keep the scripts they have written in different languages or found on the internet and put them into processes that then can be called by the workflow script.  This makes all those scripts we have written much more module.  In addition, workflow languages abstract away the submission and tracking of jobs making the program more universally transferrable between computers.

# nextflow

---

**nextflow** is a workflow framework that can be used by a bioinformatician to integrate all of her/his/their bash/python/perl/other scripts into a one cohesive pipeline that are portable, reproducible, scalable and checkpointed. Nextflow is its own DSL (Domain Specific Language) that extends a language called groovy which is an extension of Java that has language feature similarity to Python, Ruby and Smalltalk.  

  * [Introduction to nextflow](01_introductionToNextFlow.md)
  * [Creating your own workflow using nextflow](02_creatingAworkflow.md)
