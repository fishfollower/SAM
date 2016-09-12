R=R
# -> you can do    R=R-devel  make ....

PACKAGE=SAM
VERSION := $(shell sed -n '/^Version: /s///p' SAM/DESCRIPTION)

TARBALL := $(PACKAGE)_$(VERSION).tar.gz
ZIPFILE := =$(PACKAGE)_$(VERSION).zip

CPP_SRC := $(PACKAGE)/src/*.cpp

all:
	make doc-update
	make build-package
	make install
	make pdf

doc-update: $(PACKAGE)/R/*.R
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\",roclets = c(\"collate\", \"rd\"))" | $(R) --slave
	@touch doc-update

vignette-update: $(PACKAGE)/vignettes/*.Rnw
	cd $(PACKAGE)/vignettes; echo "library(knitr);knit2pdf('SAM.Rnw')" | $(R) --slave
	 mv $(PACKAGE)/vignettes/SAM.pdf $(PACKAGE)/inst/doc
	@touch vignette-update

namespace-update :: $(PACKAGE)/NAMESPACE
$(PACKAGE)/NAMESPACE: $(PACKAGE)/R/*.R
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\",roclets = c(\"namespace\"))" | $(R) --slave

build-package: $(TARBALL)
$(TARBALL): $(PACKAGE)/NAMESPACE $(CPP_SRC) $(PACKAGE)/R/*.R
	$(R) CMD build --resave-data=no $(PACKAGE)

install: $(TARBALL)
	$(R) CMD INSTALL --preclean $<
	@touch $@

quick-install: $(PACKAGE)/src/SAM.so
	$(R) CMD INSTALL $(PACKAGE)

$(PACKAGE)/src/SAM.so: $(PACKAGE)/src/SAM.cpp
	cd $(PACKAGE)/src; echo "library(TMB); compile('SAM.cpp','-O0 -g')" | $(R) --slave

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

test:
	echo "devtools::test('SAM')" | $(R) --slave

quick-check: quick-install ex-test

ex-test:
	echo "library(SAM); example(sam.fit)" | $(R) --slave

clean:
	\rm -f install doc-update $(TARBALL) $(PACKAGE).pdf $(PACKAGE)/src/SAM.so $(PACKAGE)/src/SAM.o

