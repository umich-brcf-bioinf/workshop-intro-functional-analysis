---
title: "Intro to R & RStudio: Advanced setup instructions"
author: "UM Bioinformatics Core Workshop Team"
output:
        html_document:
            theme: paper
            toc: true
            toc_depth: 6
            toc_float: true
            number_sections: false
            fig_caption: false
            markdown: GFM
            code_download: false
---

- The steps below are oriented toward advanced users who would like to install
  and configure select software on their local workstation, specifically 
  R/RStudio and dependent libraries. These steps are not
  necessary or recommended, but may provide context for advanced use-cases.

- These instructions build on the install detailed in the [basic setup
  instructions](setup_instructions.html); please see that document for more
  details on setting up Zoom, Slack, and also for getting help.

- Note that if you do not have Administrative privileges it may be difficult
  to install R/RStudio; you may need to coordinate with your
  System Admin/IT Support team to get these installed.

  #### Table of Contents
  -   [Windows setup](#windows-setup)
      -   Installing R/RStudio
      -   Notes
  -   [Macintosh setup](#macintosh-setup)
      -   Installing R/RStudio
      -   Notes
  - Installing required libraries
      - R packages
      
## Windows setup

### Installing R/RStudio (Windows)

1.1.  RStudio depends on the R programming environment, so we have to
    install that first. In a web browser, open:
    <https://cran.rstudio.com/bin/windows/base/>
    and click "Download R 4.3.2 for Windows" (the version may be
    slightly different). Open the downloaded executable to launch the R
    installer.

1.2.  The installer will walk through several options; accept all the
    defaults (by repeatedly clicking **Next**) and when prompted, click
    **Install**. The installer will show a progress bar; the process
    takes about 2 minutes; click **Finish** when prompted to close the
    installer.

1.3.  To install RStudio, in a web-browser, open:
    <https://rstudio.com/products/rstudio/download/#download>
    and click on **Download RStudio Desktop for Windows**. Open the
    downloaded executable to launch the installer.

1.4.  The installer will either prompt you to login as an Admin user or
    (if your current account has Admin privileges) simply ask you to
    allow it to make changes. Click **Yes**

1.5.  The installer will walk through several options; accept all the
    defaults (by repeatedly clicking **Next**) and when prompted, click
    **Install**. The installer will show a progress bar; the process
    takes less than one minute; click **Finish** when prompted to close
    the installer.

1.6.  Press Windows+R keys to open the **Run** dialog; type **"RStudio"**
    in the text box and hit enter. This will launch a new RStudio
    window. The RStudio window is divided into several panes. The lower
    left pane shows the **Console** tab and will show some text followed
    by a command prompt (\>):

```
R version 4.4.2 (2024-10-31 ucrt) -- "Pile of Leaves"
Copyright (C) 2024 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> 

```

1.7.  Following the workshop, you can remove R and
    RStudio. As an Admin user, go Start \> Settings \> Apps & Features.
    Click on the program to remove and click Uninstall.


## Macintosh setup

### Installing R/RStudio (Macintosh)

1.1.  RStudio depends on the R programming environment, so we have to
    install that first. In a web browser, open: 
    [https://cran.rstudio.com/bin/macosx](https://cran.rstudio.com/bin/macosx){target="_blank"}. 
    - Note that newer Macs will use Silicon build, older Macs (pre 2020) will use Intel builds. 
      (You can confirm your architecture under About this Mac: Chip; "Apple M" indicates Silicon. )
    - Right-click to download the latest version (e.g. "R-4.4.2"", though the 
    specific version may be slightly different). 
    - Open the downloaded executable to launch the R installer.

1.2.  The installer will walk through several options; accept all the
    defaults (by repeatedly clicking **Continue**) and when prompted,
    click **Install**. The installer will prompt you to confirm your
    username/password. The installer will show a progress bar; the
    process takes about 1 minutes; click **Finish** when prompted to
    close the installer.

1.3.  To install RStudio, in a web-browser, open:
-   <https://rstudio.com/products/rstudio/download/#download>
    and click on **Download RStudio Desktop for Mac**.

1.4.  Opening the downloaded executable opens a window with "RStudio" and
    your Applications folder icons. Drag the RStudio into the
    Applications folder. (If you see a dialog that claims RStudio
    already exists, click **Keep Both**.) The installer will prompt you
    to confirm your username/password.

1.5.  The installer will walk through several options; accept all the
    defaults (by repeatedly clicking **Next**) and when prompted, click
    **Install**. The installer will show a progress bar; the process
    takes less than one minute.

1.6.  When completed, open the Applications folder and double-click on the
    RStudio application. You will see a dialog
    "*RStudio.app is an app downloaded from Internet. Are you sure you
     want to open it?*"
    Click **Open**.
    This will launch a new RStudio window. The RStudio window is divided
    into several panes. The lower left pane shows the **Console** tab and
    will show some text followed by a command prompt (\>):


```
R version 4.4.2 (2024-10-31) -- "Pile of Leaves"
Copyright (C) 2024 The R Foundation for Statistical Computing
Platform: aarch64-apple-darwin20

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

    Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> 

```

1.7.  Following the workshop, you can remove R and RStudio. Open
    your Applications directory, and drag the R and RStudio application
    into the Trash.


### Installing required libraries

#### Install R packages

3.1  The workshop exercises requires the installation of special R
    packages. To install them into RStudio: 
    
    - one at a time, copy the blocks of code below and paste into the RStudio 
      **Console** tab
    - make sure the block you just pasted is highlighted and press **Enter** 
      to execute.
    - as the installation progresses, you might see red text flash by 
      and that's ok (typically an informative blurb or minor warning that has
      no downstream impact). 
      
```
required_packages = c("tidyverse",
                      "ggrepel", 
                      "RcolorBrewer")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  requireNamespace("BiocManager", quietly = TRUE)
}
BiocManager::install(required_packages, update=FALSE, ask=FALSE)

missing <- setdiff(required_packages, rownames(installed.packages()))
if (!length(missing)) {
  cat("Ready for workshop!\n")
  } else {
    cat("PROBLEM: could not install:", missing, "\n")
  }
}
```

3.5. Press **Control-Q** close RStudio; when prompted to *Save workspace
    image...*, click **Don't Save**.

#### Your workstation is ready for the workshop. Thank you for your patience and fortitude.