# Tutorials

This folder will contain any scripts, python notebooks, R markdown files, and data files for the tutorials.

**Convert ipynb to Markdown**

```{bash}
JUPYTER_NOTEBOOK=$1
jupyter nbconvert --output-dir='./docs' --to markdown ${JUPYTER_NOTEBOOK}
```

**Convert Rmd to Markdown**

Open the Rmd file in RStudio and hit the `knit` button. Move the markdown file to location in the website.
