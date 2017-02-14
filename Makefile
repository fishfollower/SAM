R=R
SHELL=/bin/bash
# -> you can do    R=R-devel  make ....

PACKAGE=stockassessment
VERSION := $(shell sed -n '/^Version: /s///p' stockassessment/DESCRIPTION)

TARBALL := $(PACKAGE)_$(VERSION).tar.gz
ZIPFILE := =$(PACKAGE)_$(VERSION).zip

CPP_SRC := $(PACKAGE)/src/*.cpp

SUBDIRS := $(wildcard tests/*/.)

.PHONY: test all $(SUBDIRS)

all:
	make doc-update
	make build-package
	make install
	make pdf

doc-update: $(PACKAGE)/R/*.R
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\")" | $(R) --slave
	@touch doc-update

vignette-update: $(PACKAGE)/vignettes/*.Rnw
	cd $(PACKAGE)/vignettes; echo "library(knitr);knit2pdf('stockassessment.Rnw')" | $(R) --slave
	mv $(PACKAGE)/vignettes/stockassessment.pdf $(PACKAGE)/inst/doc
	@touch vignette-update

namespace-update :: $(PACKAGE)/NAMESPACE
$(PACKAGE)/NAMESPACE: $(PACKAGE)/R/*.R
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\")" | $(R) --slave

build-package: $(TARBALL)
$(TARBALL): $(PACKAGE)/NAMESPACE $(CPP_SRC) $(PACKAGE)/R/*.R
	$(R) CMD build --resave-data=no $(PACKAGE)

install: $(TARBALL)
	$(R) CMD INSTALL --preclean $<
	@touch $@

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

updataData: 
	echo "library(stockassessment); \
	      source('stockassessment/tests/nscod/script.R', chdir=TRUE, echo=TRUE); \
	      nscodData <- dat; nscodConf <- conf; nscodParameters <- par; \
	      save(nscodData, file='stockassessment/data/nscodData.RData'); \
	      save(nscodConf, file='stockassessment/data/nscodConf.RData'); \
	      save(nscodParameters, file='stockassessment/data/nscodParameters.RData'); " | R --vanilla


