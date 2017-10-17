R=R
SHELL=/bin/bash

PACKAGE=stockassessment
VERSION := $(shell sed -n '/^Version: /s///p' stockassessment/DESCRIPTION)
THISSHA := $(shell git log -1 --format="%H")
TARBALL := $(PACKAGE)_$(VERSION).tar.gz
ZIPFILE := =$(PACKAGE)_$(VERSION).zip

CPP_SRC := $(PACKAGE)/src/*.cpp

.PHONY: test all updateData qi quick-install vignette-update

all:
	#make doc-update
	#make build-package
	make install
	make pdf

doc-update: $(PACKAGE)/R/*.R
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\")" | $(R) --slave
	@touch doc-update

vignette-update: vignettes/*.Rnw vignettes/*.Rmd
	cd vignettes; echo "knitr::knit2pdf('stockassessment.Rnw')" | $(R) --slave
	cd vignettes; echo "rmarkdown::render('simulate.Rmd')" | $(R) --slave
	mv vignettes/*.pdf $(PACKAGE)/inst/doc
	cd vignettes; rm -f *.{log,aux,out,tex}

namespace-update :: $(PACKAGE)/NAMESPACE
$(PACKAGE)/NAMESPACE: $(PACKAGE)/R/*.R
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\")" | $(R) --slave

build-package: $(TARBALL)
$(TARBALL): $(PACKAGE)/NAMESPACE $(CPP_SRC) $(PACKAGE)/R/*.R
	sed s/dummySHA/$(THISSHA)/g description-addon > description-addon-tmp
	mv $(PACKAGE)/DESCRIPTION old-description
	cat old-description description-addon-tmp > $(PACKAGE)/DESCRIPTION  
	rm description-addon-tmp
	$(R) CMD build --resave-data=no $(PACKAGE)
	rm $(PACKAGE)/DESCRIPTION  
	mv old-description $(PACKAGE)/DESCRIPTION  

install: $(TARBALL)
	$(R) CMD INSTALL --preclean --html $<
	@touch $@

qi:
	cd $(PACKAGE)/src; echo "library(TMB); compile('stockassessment.cpp')" | $(R) --slave
	$(R) CMD INSTALL $(PACKAGE)

quick-install: $(PACKAGE)/src/stockassessment.so
	$(R) CMD INSTALL $(PACKAGE)

$(PACKAGE)/src/stockassessment.so: $(PACKAGE)/src/stockassessment.cpp
	cd $(PACKAGE)/src; echo "library(TMB); compile('stockassessment.cpp','-O0 -g')" | $(R) --slave

unexport TEXINPUTS
pdf: $(PACKAGE).pdf
$(PACKAGE).pdf: $(PACKAGE)/man/*.Rd
	rm -f $(PACKAGE).pdf
	$(R) CMD Rd2pdf --no-preview $(PACKAGE)

build:
	$(R) CMD build $(PACKAGE)

check: $(TARBALL)
	$(R) CMD check $(TARBALL)

check-cran: $(TARBALL)
	$(R) CMD check --as-cran $(TARBALL)

quick-check: quick-install ex-test

ex-test:
	echo "library(stockassessment); example(sam.fit)" | $(R) --slave

clean:
	\rm -f install doc-update $(TARBALL) $(PACKAGE).pdf $(PACKAGE)/src/*.so $(PACKAGE)/src/*.o
	\rm -rf $(PACKAGE).Rcheck
	\rm -f stockassessment/vignettes/stockassessment.{aux,log,out,pdf,tex}

test:
	echo "devtools::test('stockassessment')" | $(R) --slave

updateData: 
	echo "library(stockassessment); \
	      source('stockassessment/tests/nscod/script.R', chdir=TRUE, echo=TRUE); \
	      nscodData <- dat; nscodConf <- conf; nscodParameters <- par; \
	      save(nscodData, file='stockassessment/data/nscodData.RData'); \
	      save(nscodConf, file='stockassessment/data/nscodConf.RData'); \
	      save(nscodParameters, file='stockassessment/data/nscodParameters.RData'); " | R --vanilla

updateDocs:
	rm -rf docs
	mkdir docs
	echo "library(Rd2md); \
	      fn<-dir('$(PACKAGE)/man'); \
	      d<-sapply(fn, function(f)Rd2markdown(paste0('$(PACKAGE)/man/',f), sub('Rd','md',paste0('docs/',f))));\
	      file.copy(paste0(find.package('$(PACKAGE)'),'/html/00Index.html'), 'docs/index.html')" | R --vanilla
	cd docs; sed -i '/<img/d' index.html
	cd docs; sed -i '/DESCRIPTION/d' index.html
	cd docs; sed -i '/User guides/d' index.html
	cd docs; sed -i 's/html/md/' index.html
	cd docs; pandoc index.html -o index.md
	cd docs; sed -i '/</d' index.md
	cd docs; sed -i '/---/d' index.md
	cd docs; sed -i 's/)/) | /' index.md
	cd docs; sed -i '/Help Pages/!{p;d;};n;a |---|---|' index.md 
	cd docs; sed -i '/Help Pages/!{p;d;};n;a | | |' index.md 
	cd docs; rm index.html

#  for later 
# 
#	-cd $(@D); echo "library(stockassessment); source('script.R'); fit.exp<-fit;save(fit.exp,file='fit.expected.Rdata') " | R --vanilla
# 
