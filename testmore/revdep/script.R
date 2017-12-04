
##' Download and build an R package from GitHub
##'
##' @param repo GitHub user and repository separated by /
##' @param ref Reference to commit or branch. Default is master
##' @param subdir Path to subdir containing the package. Should be NULL if the package is in the top directory
##' @param buildArgs Character vector of arguments passed to R CMD build. Default is "--no-build-vignettes"
##' @return Nothing
##' @author Christoffer Moesgaard Albertsen
buildFromGithub <- function(repo,
                               ref = "master",
                               subdir = NULL,
                               buildArgs = c("--no-build-vignettes")
                               ){
    topdir = tempdir()
    splitRepo <- unlist(strsplit(repo,"/"))
    if(length(splitRepo) != 2)
        stop("repo must be of the form: github_user/repository")
    user <- splitRepo[1]
    urep <- splitRepo[2]
    fil <- file.path(topdir,paste0(urep,"_",ref,".zip"))
    url <- sprintf("https://github.com/%s/archive/%s.zip",repo,ref)
    download.file(url,fil,quiet=TRUE)
    a <- utils::unzip(fil,exdir = topdir)
    descriptionPath <- a[grepl(paste0(subdir,"/DESCRIPTION"),a)]
    pkgPath <- gsub("/DESCRIPTION","",descriptionPath)
    tools::Rcmd(c("build",pkgPath,buildArgs))
}


############################################
#### Get reverse dependencies from CRAN ####
############################################
options(repos = c(CRAN = 'https://cran.r-project.org/'))   
dbUrl <- utils::available.packages(utils::contrib.url(options("repos")))
cranPkg <- unlist(tools::package_dependencies("stockassessment",
                                              db = dbUrl,
                                              which = c("Depends",
                                                        "Imports",
                                                        "LinkingTo",
                                                        "Suggests",
                                                        "Enhances"),
                                              reverse=TRUE))

#################################################################################
#### List of reverse dependencies from GitHub (arguments to buildFromGithub) ####
#################################################################################

githubPkg <- list("multiStockassessment" = list(repo="calbertsen/multi_SAM",
                                                subdir="multiStockassessment")
                  )

###############################
#### Write expected output ####
###############################


############################
#### Get package tar.gz ####
############################

if(length(cranPkg) > 0)
    sapply(cranPkg,utils::download.packages,destdir="./")

lapply(githubPkg,function(x)do.call("buildFromGithub",x))

########################
#### Check packages ####
########################


tools::check_packages_in_dir("./","--as-cran --no-build-vignettes --no-vignettes --no-manual")

checkOutput <- capture.output(tools::summarize_check_packages_in_dir_results("./"))

if(length(grep("(WARN|ERROR)",checkOutput)) > 0){
    checkRes <- "NOT_OK"
}else{
    checkRes <- "OK"
}
    



###############
#### Clean ####
###############

lf <- list.files(".",recursive=TRUE,include.dirs=TRUE)
wanted <- grepl("^(script.R|res.EXP)$",lf)
file.remove(rev(lf[!wanted]), recursive = TRUE)
        

######################
#### Write result ####
######################
cat(checkOutput,sep="\n",file="full_res.out")
cat(checkRes,"\n", file="res.out")
