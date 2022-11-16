R?=R

all: 
	@if [ -f "ignore" ]; then echo -n $(shell basename $(shell pwd)); echo "... ignored"; else time -f "%e" $(MAKE) doit; fi
doit:
	unset MAKEFLAGS; echo 'source("script.R")' | ${R} --slave > /dev/null
	@echo -n $(shell basename $(shell pwd))
	@echo -n "... "
	@diff --strip-trailing-cr res.out res.EXP; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
	  echo -n "OK "; \
	else \
	  echo "Error:"; \
	fi
	rm -f res.out
