## Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs -h` - Print help message and exit.


### You can add image as the following 

```
![Example Image](images/example_image.jpg)
```

### You can add video as the following 

```
![Example Image](images/example_video.mp4)
```

## PDF File Example

[Download the PDF](pdf/example_file.pdf)

## Embedded PDF

```
<embed src="pdf/example_file.pdf" width="800px" height="600px" />
```

### You can embed any video to mkdocs 

```
<iframe width="560" height="315" src="https://www.youtube.com/embed/dQw4w9WgXcQ" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```


### How to write equations 

1.1 Inline Equations
```
\(\sin(\omega t)\)

\\(\cos(2\omega t)\\)
```

Output of above equations is 

\(\sin(\omega t)\)

\\(\cos(2\omega t)\\)


1.2 Equation Array

```
\begin{align}
a^2 + b^2 &=& 1 \\\\
a^2 + 2\times b^2 &=& 2\times\hat{V}_1
\end{align}
```

Output is :

\begin{align}
a^2 + b^2 &=& 1 \\\\
a^2 + 2\times b^2 &=& 2\times\hat{V}_1
\end{align}


1.3 Inserting a table

```
Table-1 (contents are center-aligned)

 Method      | Description                          
 :---------: | :----------------------------------:
 `GET`       |      Fetch resource   \( \alpha \)
 `PUT`       |  Update resource
 `DELETE`    |     Delete resource

 Table-2 (contents are left-aligned)

  Method      | Description                          
  :--------- | :----------------------------------
  `GET`       |      Fetch resource  
  `PUT`       |  Update resource
  `DELETE`    |     Delete resource

  Table-3 (contents are right-aligned)

   Method      | Description                          
   ---------: | ----------------------------------:
   `GET`       |      Fetch resource  
   `PUT`       |  Update resource
   `DELETE`    |     Delete resource
```


Output is : 


Table-1 (contents are center-aligned)

 Method      | Description                          
 :---------: | :----------------------------------:
 `GET`       |      Fetch resource   \( \alpha \)
 `PUT`       |  Update resource
 `DELETE`    |     Delete resource

 Table-2 (contents are left-aligned)

  Method      | Description                          
  :--------- | :----------------------------------
  `GET`       |      Fetch resource  
  `PUT`       |  Update resource
  `DELETE`    |     Delete resource

  Table-3 (contents are right-aligned)

   Method      | Description                          
   ---------: | ----------------------------------:
   `GET`       |      Fetch resource  
   `PUT`       |  Update resource
   `DELETE`    |     Delete resource



### Font Styles 

```
Checking for \\(\bf bold\\) fonts.
```
Checking for \\(\bf bold\\) fonts.


```
Checking for \\(\it italics\\) fonts.
```
Checking for \\(\it italics\\) fonts.


```
Checking for \\(\small \text{small}\\) fonts.
```
Checking for \\(\small \text{small}\\) fonts.

```
Checking for \\(\huge \text{huge}\\) fonts.
```
Checking for \\(\huge \text{huge}\\) fonts.



```
Checking for \\( \color{red}{\text{red}}\\) fonts.
```
Checking for \\( \color{red}{\text{red}}\\) fonts.



```
Checking for \\(\huge \it \color{blue}{\text{blue}}\\) fonts.
```
Checking for \\(\huge \it \color{blue}{\text{blue}}\\) fonts.



```
Checking for \\(\Huge \tt \color{magenta}{\text{typewriter magenta}}\\) fonts.
```
Checking for \\(\Huge \tt \color{magenta}{\text{typewriter magenta}}\\) fonts.



```
Checking for \\(\Large \tt \color{green}{\mathcal{C}\mathcal{a}\mathcal{l}\mathcal{i}\mathcal{g}\mathcal{r}\mathcal{a}\mathcal{p}\mathcal{h}\mathcal{i}\mathcal{c}}\\) fonts.
```
Checking for \\(\Large \tt \color{green}{\mathcal{C}\mathcal{a}\mathcal{l}\mathcal{i}\mathcal{g}\mathcal{r}\mathcal{a}\mathcal{p}\mathcal{h}\mathcal{i}\mathcal{c}}\\) fonts.


```
Checking for \\(x_p^m \\) subscripts and superscripts.
```
Checking for \\(x_p^m \\) subscripts and superscripts.

# Example

## How to add files

This section explains how to add files to your MkDocs site.

## Step 1: Organize Your Files

Before adding files, make sure to organize them in a directory structure that makes sense for your site. For example, you might want to place images in a specific folder:


## Step 2: Add the File to Your Site

To add a file, simply place it in the appropriate directory inside your MkDocs project.

For example, to add an image:
1. Place the image in the `/assets/images/` folder.
2. Link to it in your Markdown file like this:

```markdown
![Image description](/assets/images/image.jpg)
```


## Step 3: Reference the file in your content 

```
[Download PDF](/assets/files/sample.pdf)
```



Transient Stability code

```
```