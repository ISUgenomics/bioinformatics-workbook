---
title: "Contributing"

permalink: /Contributing.html
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg

---
## How to contribute

We are really glad you're reading this! We need volunteer contributors for making the bioinformatics workbook successful! Please do not hesitate to contact us via email or issues. We want you working on the stuff that you're excited about.
{: style="text-align: justify"}

Making contribution is simple:

1. Fork the repo on GitHub
2. Clone the project to your own machine
3. Edit the files or add files using your favorite editor
4. Commit changes to your own branch
5. Push your work back up to your fork
6. Submit a Pull request so that we can review your changes

NOTE: Be sure to merge the latest from "upstream" before making a pull request!
{: style="text-align: justify"}

If you're comfortable making contributions any other way, please feel free to do it your way and send us the pull request, we will gladly review the changes.
{: style="text-align: justify"}


## How you can help

### 1. Testing/reviewing published workflows:

If you are following a particular protocol, please let us know if it worked or not. Also, we would be grateful if you let us know if changes were necessary to get it working! If you are familiar with GitHub, you can also make changes and send us the pull request. We only request that you include why the changes were necessary when submitting your pull request.
{: style="text-align: justify"}

### 2. Correcting typos or grammar:

Typos and bad grammar can make a methods section difficult to read. If you have better way to elaborate a process, we encourage your contribution. If this is the case, please make those changes and send us the pull request. Alternatively, if a section is not clear or hard to replicate you can also open a issue/bug to ask for an expansion of the existing explanation.
{: style="text-align: justify"}

### 3. Contributing to sections of published methods

In bioinformatics there are many ways to answer the same question. It is also true that using a particular method over another may influence the final results. Hence, we encourage the addition of alternative methods to the existing workflows that may be more adaptable to your organism or circumstances. Otherwise, please follow existing conventions, and explain why any methodological changes may be necessary.
{: style="text-align: justify"}

### 4. Adding new protocols/methods

Are we missing a favorite topic of yours? Well, we were just waiting for you to write that section! We are very glad that you noticed it, please add your section of interest and send us the pull request.
{: style="text-align: justify"}

### 5. Suggesting topics

If you're working on a project and you're stuck with not knowing how to proceed, then please provide us the list of topics for which you may need methods. We will add it to our list or increase the method's priority in our development of this workbook.
{: style="text-align: justify"}

### 6. Sharing your methods with us

Not time to write or add methods? No problem! Just send us your methods/protocols of doing things, either as published papers, simple writeup, weblink, wikipage, whatever the form it is. We will try to convert it and add it our workbook, crediting you! You will also be added as contributor.
{: style="text-align: justify"}

### 7. Anything missing that may make the workbook more useful?

It can be anything, please feel free to get in touch with us! We will be more than happy to hear your opinion and adjust the workbook to make it more suitable to your needs.
{: style="text-align: justify"}




## Submitting changes

Please send a [GitHub Pull Request to bioinformatics-workbook](https://github.com/ISUgenomics/bioinformatics-workbook/pull/new/master) with a clear list of what you've done (read more about [pull requests](https://help.github.com/en/articles/about-pull-requests)). Please follow our writing conventions (below) and make sure all of your commits are atomic (one feature per commit).
{: style="text-align: justify"}

Always write a clear log message for your commits. One-line messages are fine for small changes, but bigger changes should look like this:
{: style="text-align: justify"}

    $ git commit -m "A brief summary of the commit
    >
    > A paragraph describing what changed and its impact."

## Writing conventions

Start reading our code and you'll get the hang of it. We can optimize for readability.
{: style="text-align: justify"}

  * Self promotion is not allowed, unless you're _Heng Li_!
  * Use best practices whenever possible. There might be many ways to do things, but the goal here is to make it **smooth for beginners without being too confusing**.
  * To consistently find the images for each chapter, images should be saved in the `assets` folder of each section. **Please no copyrights!**
  * Methods should be generalized, but we also strongly encourage _public datasets used as examples_. Please avoid using private datasets that are not accessible, allowing others to replicate your tutorials.
  {: style="text-align: justify"}
