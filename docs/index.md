# Learning 

## This is the home page of power system lab 

For full documentation visit [mkdocs.org](https://www.mkdocs.org).

This are some example of uploaded pdf files 

<a href="/manual_sss" target="_blank">Manual SSS</a>


# YouTube Video Example

Here’s an Example embedded YouTube video:

<iframe width="560" height="315" src="https://www.youtube.com/embed/DzyX_GnSnL0?si=6_6EUCiJHiUz1dO5" frameborder="0" allowfullscreen></iframe>

# Local Video Example

Here’s an example local video embedded:

<video width="100%" height="auto" controls>
  <source src="/matrix.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>


# Screenshot Example

Following is the example of image upload

<div style="position: relative; display: inline-block;">
  <!-- Image that opens in a new tab when clicked -->
  <a href="/image.png" target="_blank">
    <img src="/image.png" alt="Screenshot" style="max-width: 30%; height: auto;">
  </a>

  <!-- Download Icon in top-right corner -->
  <a href="/image.png" download>
    <img src="https://img.icons8.com/material-outlined/24/000000/download--v1.png" 
         alt="Download Icon" 
         style="position: absolute; top: 10px; right: 10px; cursor: pointer; width: 24px;">
  </a>
</div>

## Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs -h` - Print help message and exit.

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.


## Code Annotation Examples 

### Codeblocks

Some `code` goes here

### Plain codeblock 

A plain codeblock:

```
Some code here 
def myfunction()
//some comments
```

#### Code for a specific language 

some more code with the `py` at the start:

``` py 
import tensorFlow as if 
def whatever ()
```

#### This is the title

``` py title="bubble_sort.py"
def bubble_sort(items):
    for i in range(len(items)):
    for j in range(len(items) - 1 - i):
        if items[j] > items[j + 1]:
            items[j], items[j + 1] = items[j + 1], items [j]
```

# Welcome to My Java Docs Site

This is a simple MkDocs-powered documentation site for Java tutorials.

Below is a basic Java program to print "Hello, world!" to the console:

```java
public class HelloWorld {
    public static void main(String[] args) {
        System.out.println("Hello, world!");
    }
}

Simple Tutorial: Editing and Deploying Your MkDocs Site
Step 1: Edit the Content in /docs Folder
Navigate to the /docs folder in your project directory.

Edit or add new content in Markdown files (e.g., index.md, about.md).

For example, open the index.md file and add some content like:


# Welcome to Pslab MkDocs Site
This is my awesome site built with MkDocs. 🎉

## About
Here's a brief introduction to what this site is about.
Save your changes.

Step 2: Add and Commit the Changes
After you’ve made your changes in the /docs folder, use Git to add and commit the updates.

Open your terminal or command prompt in the project directory.

Check the status to see what’s been modified:

bash
Copy
Edit
git status
Stage the changes you made to the /docs folder:

bash
Copy
Edit
git add docs/
Commit the changes with a message describing your update:

bash
Copy
Edit
git commit -m "Added new content to the site"
Step 3: Push the Changes to GitHub
Once your changes are committed, push them to your GitHub repository:

bash
Copy
Edit
git push origin main
This will push the updates to the main branch on GitHub.

Vercel will automatically detect the changes in the repo, build the site using mkdocs build, and deploy the content from the /site folder to your live website.

Step 4: Check the Live Website
After pushing to GitHub, wait a moment for Vercel to rebuild the site.

Go to your Vercel dashboard and check the deployments section for the most recent deployment.

Visit your live website and confirm that your changes have been reflected.

Recap of Required Git Commands
bash
Copy
Edit
# Check the status of your repository
git status

# Stage your changes (files you've edited)
git add docs/

# Commit your changes with a descriptive message
git commit -m "Added new content to the site"

# Push your changes to GitHub
git push origin main
Final Notes:
Always make changes in the /docs folder — this is where your content lives.

Don’t push the /site folder to GitHub! Vercel will handle that part when it builds the site.

After pushing, Vercel will automatically rebuild the site and deploy the latest version.

