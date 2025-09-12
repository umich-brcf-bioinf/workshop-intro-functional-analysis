# install.packages("remotes")
# library(remotes)
# remove.packages('klippy');
# remotes::install_github("umich-brcf-bioinf/workshop-klippy");
# install.packages("kable")
# devtools::install_github("haozhu233/kableExtra")

library(rmarkdown)
library(klippy)
setwd("~/workshop-intro-functional-analysis/")


# The html from the files below don't have the nav bar
render('source/workshop_setup/preworkshop_checklist.Rmd', output_dir='html/workshop_setup/')
render('source/workshop_setup/setup_instructions.Rmd', output_dir='html/workshop_setup/')
render('source/workshop_setup/setup_instructions_advanced.Rmd', output_dir='html/workshop_setup/')

# The html from the files below do have the nav bar, so if you make changes 
# that impact the navbar (e.g. file name changes or reordering) you should 
# re-knit all of them.

#render_site('source/instructor-cheatsheet.Rmd')
render_site('source/analysis-scripts.Rmd')
render_site('source/index.md')
render_site('source/workshop-intro.Rmd')

render_site("source/01-functional-analysis-overview.Rmd")
render_site("source/02-IntroToWebGestaltandORA.Rmd")

render_site("source/03-WebGestaltRORA.Rmd")
render_site("source/07A-advanced-visualizations.Rmd")

render_site("source/05-fcs-gsea-overview.Rmd")
render_site("source/06-WebGestaltRGSEA.Rmd")
render_site("source/04-gene-set-references.Rmd")

render_site("source/07B-advanced-visualizations.Rmd")
render_site("source/08-analysis-summary.Rmd")
render_site("source/exercises.Rmd")
render_site("source/webgestaltr-on-great-lakes.Rmd")

render_site('source/workshop-wrap-up.Rmd')

#clean_site(preview=TRUE)
