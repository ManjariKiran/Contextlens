install.packages("devtools")
install.packages("roxygen2")
devtools::create("contextlens")
devtools::document()
##https://jennybc.github.io/2014-05-12-ubc/ubc-r/session03_git.html
##Steps to connect github reporitory to local rstudio

##1)Install git
#which git
#git --version
#ruby -e "$(curl -fsSL https://raw.github.com/Homebrew/homebrew/go/install)"
#brew install git
##2) Configure git
#git config --global user.name 'ManjariKiran'
#git config --global user.email 'mkiran@github.com'
#git config --global credential.helper osxkeychain
##3)Clone git URL from github respository such as git@github.com:ManjariKiran/contextlens.git
##4) open Rstudio and create a new project
##5) choose version control
##6) Paste the address and clone the repository



