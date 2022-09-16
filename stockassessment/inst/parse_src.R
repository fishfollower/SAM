
is.namespace <- function(x){
    grepl("^[[:blank:]]*namespace ",x[1])
}

is.templated <- function(x){
    grepl("^[[:blank:]]*template[[:blank:]]*<",x[1])
}

is.struct <- function(x){
    (is.templated(x) & grepl("^[[:blank:]]*(class|struct) ",x[2])) |
        grepl("^[[:blank:]]*(class|struct) ",x[1])
}
is.enum <- function(x){
    grepl("^[[:blank:]]*enum ",x[1])
}
is.forwarddecl <- function(x)

charToLine <- function(xc, ii){
    nl <- gregexpr("\\n",xc)[[1]]
    sapply(ii, function(i) sum(nl < i)+1)
}
matchDelim <- function(x, delim, allowMismatch = FALSE, returnLine = FALSE, which.a = 1, which.b = NULL){
    xc <- paste(x,collapse="\n")
    if(!grepl(delim[1],xc) || !grepl(delim[2],xc)){
        ff <- matrix(NA,0,3)
        colnames(ff) <- c("open","close","level")
        return(ff)
    }
    a <- gregexpr(delim[1],xc)[[1]]
    b <- gregexpr(delim[2],xc)[[1]]    
    a <- a + which.a - 1
    if(is.null(which.b))
        which.b <- attr(b,"match.length")
    b <- b + which.b - 1
    v <- cbind(c(a,b),c(rep(1,length(a)),rep(-1,length(b))))
    v <- v[order(v[,1]),,drop=FALSE]
    ## Always start with open
    v0 <- min(which(v[,2]==1))
    v <- v[v0:nrow(v),,drop=FALSE]
    if(!all(cumsum(v[,2]) >= 0)){
        if(allowMismatch){
            while(!all(cumsum(v[,2]) >= 0)){
                v <- v[setdiff(seq_len(nrow(v)), min(which(cumsum(v[,2]) < 0))),,drop=FALSE]
            }
        }else{
            stop("Mismatch in delimiters")
        }
    }
    ## }
    ## Give error instead?
    ff <- matrix(NA,0,3)
    lvl <- 1
    for(i in seq_len(nrow(v))){   
        if(v[i,2]==1){                      #Open
            ff <- rbind(ff,c(open = v[i,1], close = NA, level = lvl))
            lvl <- lvl+1        
        }else{#Close
            j <- max(which(is.na(ff[,2]) & ff[,3]==lvl-1))
            ff[j,2] <- v[i,1]
            lvl <- lvl-1
        }
    }
    attr(ff,"index") <- "char"
    if(returnLine){
        ff[,1] <- charToLine(xc,ff[,1])
        ff[,2] <- charToLine(xc,ff[,2])
        attr(ff,"index") <- "line"
    }
    ff
}

clearDelim <- function(x, delim, include = FALSE, newstr = "",firstLevel = 1,skip = 0, keepEnum = FALSE, ...){
    xc <- paste(x,collapse="\n")
    ff <- matchDelim(xc, delim, ...)
    if(nrow(ff) == 0 || sum(ff[,3]>firstLevel & ff[,2]-ff[,1] > 1) == 0){
        v <- strsplit(xc,"\\n")[[1]]
        return(v[!grepl("^[[:blank:]]*$",v)])
    }
    toDelete <- ff[ff[,3]>firstLevel & # Nested more than firstLevel
                   ff[,2]-ff[,1] > 1 & # More than one line
                   ff[,1] > sum(nchar(head(x,skip))) & # After the first skip lines
                   (!grepl("^[[:blank:]]*enum[[:blank:]]+",x[charToLine(xc,ff[,1])]) | !keepEnum), # Not enum
                  ,drop=FALSE]
    
    toKeep <- rep(TRUE,nchar(xc))
    for(i in seq_len(nrow(toDelete)))
        toKeep[seq(toDelete[i,1]+(!include),toDelete[i,2]-(!include),1)] <- FALSE
    sgmnt <- cumsum(abs(c(0,diff(toKeep))))
    xspl <- strsplit(xc,"")[[1]]
    keepSgmnt <- sapply(split(toKeep,sgmnt),head,n=1)
    xL <- split(xspl,sgmnt)
    for(i in which(!keepSgmnt))
        xL[[i]] <- newstr
    xn <- do.call(paste0,lapply(xL,paste,collapse=""))
    v <- strsplit(xn,"\\n")[[1]]
    v[!grepl("^[[:blank:]]*$",v)]
}

deleteComments <- function(x){
    xn <- unlist(strsplit(x,"\\n"))
    xn2 <- gsub("//.*$","",xn)
    clearDelim(xn2, c("/\\*","\\*/"),TRUE,"",0) 
}

deleteBodies <- function(x, keepEnum = TRUE){
    xn <- unlist(strsplit(x,"\\n"))
    str <- is.struct(x)
    lvl <- as.integer(str)   
    if(str){                            # Struct or class
        xr <- clearDelim(xn, c("\\{","\\}"),FALSE,"",1,TRUE, keepEnum = TRUE)
        xr2 <- clearDelim(xr,c("[^:]:[^:]","\\{"), allowMismatch=TRUE,firstLevel = 0,which.a=2, skip = 1 + is.templated(xr), keepEnum = TRUE)
        xr <- gsub("[[:blank:]]*:?\\{\\};?",";",xr2)
    }else{                              # Functions
        xr <- clearDelim(xn, c("\\{","\\}"),TRUE,";",0, keepEnum = TRUE)
    }
    xr
}

getBody <- function(x, delim, ...){
    ff <- matchDelim(x, delim, ...)
    ii <- min(which(ff[,3]==1))
    xc <- paste(x,collapse="\n")
    unlist(strsplit(paste0(unlist(strsplit(xc,""))[seq(ff[ii,1]+1,ff[ii,2]-1,1)], collapse = ""),"\n"))
}


## classToHeader <- function(x){
##     xr <- deleteBodies(deleteComments(x))
##     xr2 <- clearDelim(xr,c(":","\\{"), allowMismatch=TRUE)
##     gsub("[[:blank:]]*:?\\{\\};?",";",xr2)
## }


## Reference classes to contain code
## File (List of namespaces & CodeChunks)
## CodeChunk (members: file, original lines, src code; functions: toHeader, Specialize, SpecializeHeader, Instantiate)
CodeChunk <- setRefClass("CodeChunk",
                         fields = list(file = "character",
                                       fileLines = "integer",
                                       src = "character"),
                         methods = list(
                             show = function(){
                                 cat("From file:",file,"\n")
                                 nd <- ceiling(log10(fileLines[2]))
                                 ln <- sprintf(sprintf("%%0%dd",nd),as.integer(seq(fileLines[1],fileLines[2],1)))
                                 cat(paste0(ln,": ",src),sep="\n")
                             },
                             writeCpp = function(newfile = "", append = TRUE){
                                 cat(src,sep="\n",file=newfile, append = append)
                             },
                             writeHpp = function(newfile = "", append = TRUE){
                                 srcH = deleteBodies(src, TRUE)
                                 cat(srcH,sep="\n",file=newfile, append = append)
                             },
                             templated = function(){ is.templated(src) },
                             getIndex = function() { fileLines[1] }
                         )
                         )

## - Namespace (List of CodeChunk)
CodeChunk_Namespace <- setRefClass("CodeChunk_Namespace",
                                   fields = list(objects = "list"),
                                   contains = "CodeChunk",
                                   methods = list(
                                       show = function(){
                                           cat("From file:",file,"\n")
                                           nd <- ceiling(log10(fileLines[2]))
                                           ln <- sprintf(sprintf("%%0%dd",nd),as.integer(c(fileLines[1],fileLines[2])))
                                           cat(paste0(ln[1]," - ",ln[2],": ",src),sep="\n")
                                           for(i in seq_along(objects))
                                               objects[[i]]$show()
                                       },
                                       writeCpp = function(newfile = "", append = TRUE){
                                           cat(src,"{",sep="\n",file=newfile, append = append)
                                           for(i in seq_along(objects))
                                               objects[[i]]$writeCpp(newfile,append)
                                           cat(src,"}",sep="\n",file=newfile, append = append)
                                       },
                                       writeHpp = function(newfile = "", append = TRUE){
                                           cat(src,"{",sep="\n",file=newfile, append = append)
                                           for(i in seq_along(objects))
                                               objects[[i]]$writeHpp(newfile,append)
                                           cat(src,"}",sep="\n",file=newfile, append = append)
                                       })
                                   )
CodeChunk_NamespaceOmit <- setRefClass("CodeChunk_NamespaceOmit",
                                   contains = "CodeChunk",
                                   methods = list(                                     
                                       writeHpp = function(newfile = "", append = TRUE){
                                           ## cat(src,"{",sep="\n",file=newfile, append = append)
                                           ## for(i in seq_along(objects))
                                           ##     objects[[i]]$writeHpp(newfile,append)
                                           ## cat(src,"}",sep="\n",file=newfile, append = append)
                                       })
                                   )
                                   

## - Preprocessor (should prbably check for Type?)
CodeChunk_Preprocessor <- setRefClass("CodeChunk_Preprocessor",
                                      fields = list(isHeaderGuard = "logical"),
                                      contains = "CodeChunk",
                                      methods = list(
                                          writeCpp = function(newfile = "", append = TRUE){
                                              if(!isHeaderGuard)
                                                  cat(src,sep="\n",file=newfile, append = append)
                                          },
                                          writeHpp = function(newfile = "", append = TRUE){
                                              if(!isHeaderGuard)
                                                  cat(src,sep="\n",file=newfile, append = append)
                                          },
                                          templated = function() { FALSE }
                                      )
                                      )

## - Struct
CodeChunk_Struct <- setRefClass("CodeChunk_Struct",
                                contains = "CodeChunk"
                                ## To Cpp, this should only output member functions:
                                ## template<class Type, T>
                                ## T Struct<Type>::functionName(T a){...};
                         )
## - Function
CodeChunk_Function <- setRefClass("CodeChunk_Function",
                                contains = "CodeChunk"
                         )                        
## - Enum
CodeChunk_Enum <- setRefClass("CodeChunk_Enum",
                                contains = "CodeChunk",
                                      methods = list(
                                          writeCpp = function(newfile = "", append = TRUE){
                                              ##cat(src,sep="\n",file=newfile, append = append)
                                          },
                                          writeHpp = function(newfile = "", append = TRUE){
                                              cat(src,sep="\n",file=newfile, append = append)
                                          },
                                          templated = function() { FALSE }
                                      )
                         )                        



splitCode <- function(x, file = NA_character_, recursive = FALSE){
    x <- x[!grepl("^[[:blank:]]*$",x)]
    iRemain <- seq_along(x)
    xr <- x[iRemain]
    ## Find namespaces
    nni <- grep("^[[:blank:]]*namespace[[:blank:]]",xr)
    xNm <- list()
    if(length(nni) > 0){                # There are namespaces
        ## Find namespace delimiter
        ff <- matchDelim(xr,c("\\{","\\}"),returnLine=TRUE)
        delimOnLine <- grepl("\\{[[:space:]]*$",xr[nni])
        ff <- ff[match(nni - ifelse(delimOnLine,0,1),ff[,1]),,drop=FALSE]
        ## Check that they are in first level
        ff <- ff[ff[,3]==1,,drop=FALSE]
        if(nrow(ff) > 0){               # There are namespaces to handle
            xNm <- lapply(seq_len(nrow(ff)), function(i){
                nscode <- xr[nni[i]:ff[i,2]]
                if(recursive){
                    CodeChunk_Namespace$new(file = as.character(file),
                                            fileLines = as.integer(c(ff[i,1],ff[i,2])),
                                            src = gsub("[[:blank:]]*\\{[[:blank:]]*","",xr[ff[i,1]]),
                                            objects = splitCode(getBody(nscode,c("\\{","\\}"))))
                }else{
                    CodeChunk_Namespace$new(file = as.character(file),
                                            fileLines = as.integer(c(ff[i,1],ff[i,2])),
                                            src = nscode)
                }
            })
            iRemain <- setdiff(iRemain,iRemain[unlist(lapply(xNm,function(cc){ seq(cc$fileLines[1],cc$fileLines[2],1) }))])
            xr <- x[iRemain]
        } 
    }
    ## Find preprocessor flags (and remove before next step)
    ppi <- grep("^[[:blank:]]*#.+[^\\\\]?$",xr)
    ppc <- grep("\\\\$",xr)
    isGuard <- rep(FALSE,length(ppi))
    if(grepl("#pragma once",xr[ppi[1]]))
        isGuard[1] <- TRUE
    if(grepl("#ifndef ",xr[ppi[1+isGuard[1]]]) &&
       grepl("#define ",xr[ppi[1+isGuard[1]+1]]) &&
       grepl("#endif",xr[tail(ppi,1)])){
        isGuard[c(1+isGuard[1],1+isGuard[1]+1,length(ppi))] <- TRUE
    }
    xPm <- lapply(ppi, function(i){
        if(grepl("\\\\$",x[i])){         #Multi line
            if(length(ppc[ppc>=i]) == 1){
                to <- i+1
            }else{
                j <- min(which(diff(ppc[ppc>=i])>1))-1
                to <- ppc[ppc>=i][j+1]+1
            }
        }else{
            to <- i
        }
        CodeChunk_Preprocessor$new(file = as.character(file), fileLines = as.integer(c(i,to)), src = x[seq(i,to,1)], isHeaderGuard = isGuard[match(i,ppi)])
    })

    ## Find 
    
    ## Remove preprocessor code
    iRemain <- setdiff(iRemain,iRemain[unlist(lapply(xPm,function(cc){ seq(cc$fileLines[1],cc$fileLines[2],1) }))])
    xr <- x[iRemain]
    
    ## Find structs / functions / namespaces
    ff <- matchDelim(xr,c("\\{","\\}"),FALSE,TRUE)
    ffIn <- ff[ff[,3]==1,,drop=FALSE]
    xL <- lapply(seq_len(nrow(ffIn)),function(i){
        from <- ffIn[i,1]
        if(grepl("^[[:blank:]]*template[[:blank:]]*<",xr[pmax(1,from-1)]))
            from <- from-1
        to <- ffIn[i,2]
        r <- xr[seq(from,to,1)]
        if(is.struct(r)){
            return(CodeChunk_Struct$new(file = file, fileLines = iRemain[c(from,to)], src = r))
        }else if(is.enum(r)){
            return(CodeChunk_Enum$new(file = file, fileLines = iRemain[c(from,to)], src = r))
        }
        return(CodeChunk_Function$new(file = file, fileLines = iRemain[c(from,to)], src = r))
    })
    xOut <- c(xNm,xPm,xL)
    sl <- sapply(xOut, function(xx) xx$getIndex())
    xOut[order(sl)]    
}

parseSAMFile <- function(f){
    cat("Parsing",f,"\n")
    l <- readLines(file.path("inst","include",f))
    l <- deleteComments(l)
    l <- l[!grepl("^[[:blank:]]*$",l)]
    xx <- splitCode(l,f)
    fC <- file.path("src",paste0("A_",gsub("\\.hpp",".cpp",f)))
    file.create(fC)
    fH <- file.path("src",paste0("A_",gsub("\\.hpp",".h",f)))
    file.create(fH)
    ## Make header guards and includes
    cat(sapply(c("#pragma once","#ifndef _SAM_AUTOGENERATED_%s_","#define _SAM_AUTOGENERATED_%s_"),sprintf,toupper(gsub("\\.hpp","",f))),file=fH,sep="\n",append=TRUE)
    cat(c("#include \"TMB.h\"","#include \"A_define.h\""),file=fH,sep="\n",append=TRUE)
    ## Include in cpp
    cat(sprintf("#include \"%s\"",paste0("A_",gsub("\\.hpp",".h",f))),file=fC,sep="\n",append=TRUE)
    for(i in seq_along(xx)){
        xx[[i]]$writeCpp(fC)
        xx[[i]]$writeHpp(fH)
    }
    cat("#endif",file=fH,sep="\n",append=TRUE)
}

makeSAMinclude <- function(){
    f <- list.files("src","A_.+\\.h")
    cat(paste("#include",f),file="src/SAM.h",sep="\n")
}
