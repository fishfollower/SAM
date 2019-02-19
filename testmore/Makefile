R?=R

all: 
	@if [ -f "ignore" ]; then echo -n $(shell basename $(shell pwd)); echo "... ignored"; else $(MAKE) doit; fi
doit:
	unset MAKEFLAGS; echo 'source("script.R")' | ${R} --slave > /dev/null
	@echo -n $(shell basename $(shell pwd))
	@echo -n "... "
	@diff --strip-trailing-cr res.out res.EXP; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
	  echo "OK"; \
	else \
	  echo "Error:"; \
	fi
	rm -f res.out
