## Dynamic websites

Dynamic websites pull information from a database to fill in the content on a webpage. for example, the search results page you are shown didn’t already exist as a full HTML page; instead, Amazon.com has a template for search results page that includes things all results pages share (like the main menu and Amazon logo), but it queries the database to insert the results of that search you initiated into that template.

Static websites, on the other hand, do not use a database to store information; instead, all information to be displayed on each webpage is already contained in an HTML file for that webpage. The HTML pages that make up a static site can be completely written by hand, or you can offload some of this work using something like Jekyll.

Jekyll is software that helps you “generate” or create a static website (you may see Jekyll described as a “static site generator”). Jekyll takes page templates—those things like main menus and footers that you’d like shared across all the web pages on your site, where manually writing the HTML to include them on every webpage would be time-consuming. These templates are combined with other files with specific information (e.g. a file for each blog post on the site) to generate full HTML pages for website visitors to see. Jekyll doesn’t need to do anything like querying a database and creating a new HTML page (or filling in a partial one) when you visit a webpage; it’s already got the HTML pages fully formed, and it just updates them when/if they ever change.

### GitHub & GitHub Pages

GitHub Pages is a free place to store the files that run a website and host that website for people to visit (it only works for particular types of website, like basic HTML sites or Jekyll sites, and does not host databases).

GitHub is a visual way to use git, a system for versioning: keeping track of changes to computer files (including code and text documents) over time (as explained above). If you’re curious, here’s a friendly lesson for exploring GitHub.

##### Built for blogging:
 Jekyll was built to support blog posts, so it’s easy to blog (add new, date-sorted content) and do related tasks like display an archive of all blog posts by month, or include a link to the three most recent blog posts at the bottom of each post.

## _ site folder
This is where the generated site will be placed (by default) once Jekyll is done transforming it. It’s probably a good idea to add this to your .gitignore file.
