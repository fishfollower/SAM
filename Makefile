R?=R

PACKAGE=stockassessment
VERSION := $(shell sed -n '/^Version: /s///p' stockassessment/DESCRIPTION)
THISSHA := $(shell git log -1 --format="%H")
TARBALL := $(PACKAGE)_$(VERSION).tar.gz
ZIPFILE := =$(PACKAGE)_$(VERSION).zip

CPP_SRC := $(PACKAGE)/src/*.cpp $(PACKAGE)/inst/include/*.hpp

SUBDIRS := $(wildcard testmore/*/.)

ifeq (testonemore,$(firstword $(MAKECMDGOALS)))
  ARG := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(ARG):;@:)
endif

ifeq (webtestone,$(firstword $(MAKECMDGOALS)))
  ARG := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(ARG):;@:)
endif

ifeq (webtest,$(firstword $(MAKECMDGOALS)))
  ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(ARGS):;@:)
endif

testfiles := $(foreach dir,$(ARGS),$(dir)/OK)

.PHONY: webtestfromfile webtestone webtest test testmoreseq testonemore testmore $(SUBDIRS) all updateData qi quick-install vignette-update 

all:
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

$(PACKAGE)/src/stockassessment.so: $(PACKAGE)/src/stockassessment.cpp $(CPP_SRC)
	touch $(PACKAGE)/src/stockassessment.cpp
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
	rm -f install doc-update $(TARBALL) $(PACKAGE).pdf $(PACKAGE)/src/*.so $(PACKAGE)/src/*.o
	rm -rf $(PACKAGE).Rcheck
	rm -f stockassessment/vignettes/stockassessment.{aux,log,out,pdf,tex}

test:
	echo "devtools::test('stockassessment')" | $(R) --slave

updateData: 
	echo "library(stockassessment); \
	      source('stockassessment/tests/nscod/script.R', chdir=TRUE, echo=TRUE); \
	      nscodData <- dat; nscodConf <- conf; nscodParameters <- par; \
	      save(nscodData, file='stockassessment/data/nscodData.RData'); \
	      save(nscodConf, file='stockassessment/data/nscodConf.RData'); \
	      save(nscodParameters, file='stockassessment/data/nscodParameters.RData'); " | R --slave

updateDocs:
	rm -rf docs
	mkdir docs
	echo "library(Rd2md); \
	      fn<-dir('$(PACKAGE)/man'); \
	      d<-sapply(fn, function(f)Rd2markdown(paste0('$(PACKAGE)/man/',f), sub('Rd','md',paste0('docs/',f))));\
	      file.copy(paste0(find.package('$(PACKAGE)'),'/html/00Index.html'), 'docs/index.html')" | R --slave
	R CMD Rdconv -t html $(PACKAGE)/man/sam.fit.Rd -o temp.html
	sed -i '/page for sam.fit/d' temp.html
	pandoc temp.html -t markdown_github -o docs/sam.fit.md; rm temp.html;
	cd docs; sed -i '/<img/d; /DESCRIPTION/d; /User guides/d; s/html/md/' index.html
	cd docs; pandoc index.html -t markdown_github -o index.md
	cd docs; rm index.html

NPROCS:=1
OS:=$(shell uname -s)

ifeq ($(OS),Linux)
  NPROCS:=$(shell grep -c ^processor /proc/cpuinfo)
endif

MAKEFLAGS += --silent

testmore:
	$(MAKE) -j $(NPROCS) testmoreseq

testmoreseq: $(SUBDIRS)

testonemore:
	@$(MAKE) testmore/$(ARG)/.

$(SUBDIRS):
	@cp testmore/Makefile $@
	@$(MAKE) -i -s -C $@
	@rm -f $@/Makefile

webtestfromfile:
	$(MAKE) -j $(NPROCS) webtest $$(< webteststocklist)

webtest: $(testfiles)

webtestone:
	@wget -q -r -np -nH --cut-dirs=4 -R index.html* https://www.stockassessment.org/datadisk/stockassessment/userdirs/user3/$(ARG)/
	@sed -i 's/useR = Rnewest/useR = R/' $(ARG)/Makefile
	@sed -i 's/--vanilla/ /' $(ARG)/Makefile
	@mv $(ARG)/run/model.RData $(ARG); 
	@rm -f $(ARG)/run/curver
	@touch $(ARG)/data/*
	@$(MAKE) -s -C $(ARG) model
	@echo "load('$(ARG)/model.RData'); old<-fit[['pl']]; \
	       load('$(ARG)/run/model.RData'); new<-fit[['pl']];\
	       cat('$(ARG)...',ifelse(all.equal(old,new,check.attributes=FALSE),'OK','FAIL'),'\n')"   | R --slave
	@touch $(ARG)/OK

$(testfiles):
	@$(MAKE) -s webtestone $(@D)
	@rm -rf $(@D)
