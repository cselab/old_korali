# Dependencies
Use pip install
1. mkdocs
2. [mkdocs-material](https://squidfunk.github.io/mkdocs-material/getting-started/)
3. pymdown-extensions

# mkdocs

1. mkdocs gh-deploy : make the site and push to the repo
2. mkdocs serve : run a local server to http://127.0.0.1:8000


# images:

If you want to scale an image you have to include directly html code. For example,
```
<img src="../images/psi.pdf" width="80%" height="80%"/>
```
The folder `images` is located in the same folder with the file the md file. The `../` is needed for the html final site.
