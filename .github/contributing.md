# How to contribute

We are really glad you're reading this! We need volunteer contributors for making this endeavour successful! Please don't hesitate to contact us via email or issues, 
we want you working on the stuff that you're excited about.

To begin, here are links to answer some of the common questions:

  * Roadmap (coming soon)
  * Our vision (coming soon)
  * Funding agencies (coming soon)
  * Issues and bugs

# How can you help:

1. Testing the published workflows:

If you are following a particular protocol, please let us know if it worked or not. Also, we would be grateful if you let us know 
what changes you had to make to get it working! If you are familiar with GitHub, you also make those changes and send us the pull request.
We only request that you include why you had to make that changes for it to work in your pull request.

2. Correcting typos or grammar:
We know that typos and bad grammar makes the method section very confusing. If you have better words to explain exiiting methods, that otherwise
is very confusing, please make those changes and send us the pull request. Alternatively, if a section is not understandable or hard to replicate
you can also open a issue/bug asking us for more explaination.

3. Contributing to sections of the published methods

We understand that in Bioinormatics there is more than one way to do things. It is especially true selcting a particular method might
infulence the final outcome (results). Hence, we encourage you to add alternative methods for the exisiting workflows that you think might be 
more suitable (for particular organism or under particular circumstances). Please follow the exisiting conventions and explain it why they
should be using other mehtods instead.

4. Adding new protocol/methods

Are we missing the favorite topic of yours? Well, we were just waiting for you to write that section becuse we didn't wanted to have
any mistakes in your favorite topic! We are very glad that you noticed it, please go ahead and add the section and send us the pull request, we will
gladly accept. 

5. Suggesting topics

If you're working on a project and you're stuck with not knowing how to proceed, then please provide us the list of topics that you need 
methods for. We will add it to our list and try to get that up for you to test as soon as our time permits.

6. Anything that is not covered here, but will make the workbook more useful!






## Testing

We have a handful of Cucumber features, but most of our testbed consists of RSpec examples. Please write RSpec examples for new code you create.

## Submitting changes

Please send a [GitHub Pull Request to opengovernment](https://github.com/opengovernment/opengovernment/pull/new/master) with a clear list of what you've done (read more about [pull requests](http://help.github.com/pull-requests/)). When you send a pull request, we will love you forever if you include RSpec examples. We can always use more test coverage. Please follow our coding conventions (below) and make sure all of your commits are atomic (one feature per commit).

Always write a clear log message for your commits. One-line messages are fine for small changes, but bigger changes should look like this:

    $ git commit -m "A brief summary of the commit
    > 
    > A paragraph describing what changed and its impact."

## Coding conventions

Start reading our code and you'll get the hang of it. We optimize for readability:

  * We indent using two spaces (soft tabs)
  * We use HAML for all views
  * We avoid logic in views, putting HTML generators into helpers
  * We ALWAYS put spaces after list items and method parameters (`[1, 2, 3]`, not `[1,2,3]`), around operators (`x += 1`, not `x+=1`), and around hash arrows.
  * This is open source software. Consider the people who will read your code, and make it look nice for them. It's sort of like driving a car: Perhaps you love doing donuts when you're alone, but with passengers the goal is to make the ride as smooth as possible.
  * So that we can consistently serve images from the CDN, always use image_path or image_tag when referring to images. Never prepend "/images/" when using image_path or image_tag.
  * Also for the CDN, always use cwd-relative paths rather than root-relative paths in image URLs in any CSS. So instead of url('/images/blah.gif'), use url('../images/blah.gif').

Thanks,
Carl Tashian, Participatory Politics Foundation
