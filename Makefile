R?=R

PACKAGE=stockassessment
VERSION := $(shell sed -n '/^Version: /s///p' stockassessment/DESCRIPTION)
THISSHA := $(shell git log -1 --format="%H")
TARBALL := $(PACKAGE)_$(VERSION).tar.gz
ZIPFILE := $(PACKAGE)_$(VERSION).zip

CPP_SRC := $(PACKAGE)/src/*.cpp $(PACKAGE)/src/*.h $(PACKAGE)/inst/include/*.hpp $(PACKAGE)/inst/include/SAM/*.hpp
R_FILES := $(PACKAGE)/R/*.R

SUBDIRS := $(wildcard testmore/*/.)

ifeq (testonemore,$(firstword $(MAKECMDGOALS)))
  ARG := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(ARG):;@:)
endif

ifeq (testoneleak,$(firstword $(MAKECMDGOALS)))
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

NPROCS:=1
OS:=$(shell uname -s)

ifeq ($(OS),Linux)
  #NPROCS:=$(shell grep -c ^processor /proc/cpuinfo)
# Use number of physical cores -1 instead of number of threads
  NPROCS:=$(shell grep 'cpu cores' /proc/cpuinfo | uniq | sed "s/.*: //g" | xargs -n 1 expr -1 + )
endif



.PHONY: webtestfromfile webtestone webtest test testmoreseq testonemore testmore $(SUBDIRS) all updateData qi quick-install vignette-update 

all:
	make install

$(PACKAGE)/configure:  $(PACKAGE)/configure.ac
	cd $(PACKAGE) && autoconf

doc-update: $(R_FILES)
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\")" | $(R) --slave
	sed -i /RoxygenNote/d $(PACKAGE)/DESCRIPTION
	@touch doc-update

vignette-update: vignettes/*.Rnw vignettes/*.Rmd
	cd vignettes; echo "knitr::knit2pdf('stockassessment.Rnw')" | $(R) --slave
	cd vignettes; echo "rmarkdown::render('simulate.Rmd')" | $(R) --slave
	mv vignettes/*.pdf $(PACKAGE)/inst/doc
	cd vignettes; rm -f *.{log,aux,out,tex}

namespace-update :: $(PACKAGE)/NAMESPACE
$(PACKAGE)/NAMESPACE: $(R_FILES)
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\")" | $(R) --slave

build-package: $(TARBALL)
$(TARBALL): $(PACKAGE)/NAMESPACE $(CPP_SRC) $(PACKAGE)/R/*.R $(PACKAGE)/configure
	sed s/dummySHA/$(THISSHA)/g description-addon > description-addon-tmp
	mv $(PACKAGE)/DESCRIPTION old-description
	sed -i -E '/^(Remote|Github)/d' old-description
	cat old-description description-addon-tmp > $(PACKAGE)/DESCRIPTION  
	rm description-addon-tmp
	$(R) CMD build --no-manual --resave-data=no $(PACKAGE)
	rm $(PACKAGE)/DESCRIPTION  
	mv old-description $(PACKAGE)/DESCRIPTION  
	sed -i /RoxygenNote/d $(PACKAGE)/DESCRIPTION

install: $(TARBALL)
	MAKE="$(MAKE) -j$(NPROCS)" $(R) CMD INSTALL --preclean --html $<
	@touch $@

qi:	$(PACKAGE)/configure
	$(error This will not work as intended. Use 'make quick-install' instead)
	 cd $(PACKAGE)/src; echo "library(TMB); compile('stockassessment.cpp')" | $(R) --slave
	MAKE="$(MAKE) -j$(NPROCS)" $(R) CMD INSTALL $(PACKAGE)

quick-install: $(CPP_SRC) $(PACKAGE)/configure $(R_FILES)
	$(info Running R CMD INSTALL with ${NPROCS} threads)
	MAKE="$(MAKE) -j$(NPROCS)" $(R) CMD INSTALL $(PACKAGE)

quick-install-debug: $(CPP_SRC) $(PACKAGE)/configure $(R_FILES)
	$(R) CMD INSTALL --configure-args='--enable-debug' $(PACKAGE)

quick-install-analyzer: $(CPP_SRC) $(PACKAGE)/configure $(R_FILES)
	$(R) CMD INSTALL --configure-args='--enable-debug --enable-analyze' $(PACKAGE)


$(PACKAGE)/src/stockassessment.so: $(PACKAGE)/src/stockassessment.cpp $(CPP_SRC) $(PACKAGE)/configure
	$(error This will not work as intended. Use 'make quick-install' instead)
	touch $(PACKAGE)/src/stockassessment.cpp
	cd $(PACKAGE)/src; echo "library(TMB); compile('stockassessment.cpp','-O0 -g', libinit=FALSE)" | $(R) --slave

unexport TEXINPUTS
pdf: $(PACKAGE).pdf
$(PACKAGE).pdf: $(PACKAGE)/man/*.Rd
	rm -f $(PACKAGE).pdf
	$(R) CMD Rd2pdf --no-preview $(PACKAGE)

build:
	MAKE="$(MAKE) -j$(NPROCS)" $(R) CMD build $(PACKAGE)

check: $(TARBALL)
	MAKE="$(MAKE) -j$(NPROCS)" $(R) CMD check $(TARBALL)

check-cran: $(TARBALL)
	MAKE="$(MAKE) -j$(NPROCS)" $(R) CMD check --as-cran $(TARBALL)

quick-check: quick-install ex-test

ex-test:
	echo "library(stockassessment); example(sam.fit)" | $(R) --slave

clean:
	cd $(PACKAGE) && ./cleanup && cd ..
	rm -f install doc-update $(TARBALL) $(PACKAGE).pdf $(PACKAGE)/src/*.so $(PACKAGE)/src/*.o
	rm -rf $(PACKAGE).Rcheck
	rm -f stockassessment/vignettes/stockassessment.{aux,log,out,pdf,tex}

test:
	echo "devtools::test('stockassessment')" | $(R) --slave

updateData: 
	echo "library(stockassessment); \
	      source('stockassessment/tests/nscod/script.R', chdir=TRUE, echo=TRUE); \
	      attr(par,'what') <- NULL; par$missing <- NULL; \
              nscodData <- dat; nscodConf <- conf; nscodParameters <- par; \
	      save(nscodData, file='stockassessment/data/nscodData.RData', version=2); \
	      save(nscodConf, file='stockassessment/data/nscodConf.RData', version=2); \
	      save(nscodParameters, file='stockassessment/data/nscodParameters.RData', version=2); " | $(R) --slave

updateDocs:
	rm -rf docs
	mkdir docs
	echo "library(Rd2md); \
	      fn<-dir('$(PACKAGE)/man'); \
	      d<-sapply(fn, function(f)Rd2markdown(paste0('$(PACKAGE)/man/',f), sub('Rd','md',paste0('docs/',f))));\
	      file.copy(paste0(find.package('$(PACKAGE)'),'/html/00Index.html'), 'docs/index.html')" | $(R) --slave
	$(R) CMD Rdconv -t html $(PACKAGE)/man/sam.fit.Rd -o temp.html
	sed -i '/page for sam.fit/d' temp.html
	pandoc temp.html -t markdown_github -o docs/sam.fit.md; rm temp.html;
	cd docs; sed -i '/<img/d; /DESCRIPTION/d; /User guides/d; s/html/md/' index.html
	cd docs; pandoc index.html -t markdown_github -o index.md
	cd docs; rm index.html


MAKEFLAGS += --silent

testmore:
	$(MAKE) -j $(NPROCS) testmoreseq

testmoreseq: $(SUBDIRS)

testonemore:
	@$(MAKE) testmore/$(ARG)/.

testoneleak:
	$(R) -d 'valgrind --leak-check=full --show-leak-kinds=all' -f testmore/$(ARG)/script.R > valgrind_${ARG}.txt

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
	@echo "load('$(ARG)/model.RData'); old<-tail(summary(fit),1); \
	       load('$(ARG)/run/model.RData'); new<-tail(summary(fit),1);\
	       cat('$(ARG)...',ifelse(all.equal(old,new,check.attributes=FALSE),'OK','FAIL'),'\n')"   | R --slave
	@touch $(ARG)/OK

$(testfiles):
	@$(MAKE) -s webtestone $(@D)
	@rm -rf $(@D)

