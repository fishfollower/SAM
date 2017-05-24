xtab<-function(x,caption='Table X.', file=stdout(), width='"100%"', cornername='', dec=rep(1,ncol(x))){
  nc<-ncol(x)
  lin<-paste('<table width=',width,'>', sep='')
  lin<-c(lin,sub('$','</td></tr>',sub('\\. |\\.$','.</b> ',
         sub('^', paste('<tr><td colspan=',nc+1,'><b>',sep=''), caption))))
  hr<-paste('<tr><td colspan=',nc+1,'><hr noshade></td></tr>', sep='')
  lin<-c(lin,hr)
  cnames<-colnames(x)
  cnames<-paste(sub('$','</b></td>',sub('^','<td align=right><b>',cnames)), collapse='\t')
  lin<-c(lin,paste('<tr>',paste('<td align=left><b>',cornername,'</b></td>',sep=''),cnames,'</tr>'))
  lin<-c(lin,hr)
  rnames<-sub('$','</b></td>',sub('^','<tr> <td align=left><b>',rownames(x)))
  xx<-sapply(1:ncol(x),function(i)sub('NA','  ',formatC(round(x[,i,drop=FALSE],dec[i]),digits=dec[i], format='f')))
  x<-relist(xx, skeleton=x)
  for(i in 1:nrow(x)){
    thisline<-paste(rnames[i],paste(sub('$','</td>',sub('^','<td align=right>',x[i,])), collapse='\t'),'</tr>', sep='')
    lin<-c(lin,thisline)
  }
  lin<-c(lin,hr)
  lin<-c(lin,'</table><br>\n')
  writeLines(lin,con=file)
}

xtab2<-function(x,caption='Table X.', file=stdout(), width='"100%"', cornername='', dec=matrix(1,ncol=ncol(x),nrow=nrow(x))){
  nc<-ncol(x)
  lin<-paste('<table width=',width,'>', sep='')
  lin<-c(lin,sub('$','</td></tr>',sub('\\. |\\.$','.</b> ',
         sub('^', paste('<tr><td colspan=',nc+1,'><b>',sep=''), caption))))
  hr<-paste('<tr><td colspan=',nc+1,'><hr noshade></td></tr>', sep='')
  lin<-c(lin,hr)
  cnames<-colnames(x)
  cnames<-paste(sub('$','</b></td>',sub('^','<td align=right><b>',cnames)), collapse='\t')
  lin<-c(lin,paste('<tr>',paste('<td align=left><b>',cornername,'</b></td>',sep=''),cnames,'</tr>'))
  lin<-c(lin,hr)
  rnames<-sub('$','</b></td>',sub('^','<tr> <td align=left><b>',rownames(x)))
  xx<-x
  for(i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      x[i,j]<-sub('NA','  ',formatC(round(xx[i,j],dec[i,j]),digits=dec[i,j], format='f'))
    }
  }

  for(i in 1:nrow(x)){
    thisline<-paste(rnames[i],paste(sub('$','</td>',sub('^','<td align=right>',x[i,])), collapse='\t'),'</tr>', sep='')
    lin<-c(lin,thisline)
  }
  lin<-c(lin,hr)
  lin<-c(lin,'</table><br>\n')
  writeLines(lin,con=file)
}

dummyplot<-function(text='This plot is intentionally left blank'){
  plot(c(0,1),c(0,1),axes=FALSE, xlab='', ylab='' , type='n')
  box()
  text(.5,.5,labels=text)
}

stampit<-function() {
  x<-system('svn info', intern=TRUE)
  udir<-sub('/datadisk/stockassessment/userdirs/','',getwd())
  udir<-sub('user[[:digit:]]+/','',udir)
  udir<-sub('/res','',udir)
  txt<-paste("stockassessment.org",udir,sub("Revision: ", "r",x[grep("Revision",x)]), sep=', ')
  ## Function modified from Frank Harrell's Hmisc library 
  stamp <- function(string = "", print = TRUE, plot = TRUE) {
    opar <- par(yaxt = "s", xaxt = "s", xpd = NA)
    on.exit(par(opar))
    plt <- par("plt")
    usr <- par("usr")
    xcoord <- usr[2] + (usr[2] - usr[1])/(plt[2] - plt[1]) * (1 - plt[2]) - 0.6 * strwidth("m")
    ycoord <- usr[3] - diff(usr[3:4])/diff(plt[3:4]) * (plt[3]) + 0.6 * strheight("m")
    if (par("xlog")) xcoord <- 10^(xcoord)
    if (par("ylog")) ycoord <- 10^(ycoord)
    text(xcoord, ycoord, string, adj = 1)
    invisible(string)
  }

  oldpar <- par(mfrow = c(1, 1), cex = 0.5)
  on.exit(par(oldpar))
  stamp(string = txt, print = FALSE, plot = TRUE)
  invisible()
}

plotcounter<-1 
tit.list<-list()
setcap<-function(title="", caption=""){   
 tit.list[length(tit.list)+1]<<-paste("# Title",plotcounter)
 tit.list[length(tit.list)+1]<<-paste(title)
 tit.list[length(tit.list)+1]<<-paste("# Caption",plotcounter)
 tit.list[length(tit.list)+1]<<-paste(caption)
 plotcounter<<-plotcounter+1 
}
