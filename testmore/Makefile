all: 
	@echo 'source("script.R")' | R --slave > /dev/null
	@echo -n $(shell basename $(shell pwd))
	@echo -n "... "
	@cmp -s res.out res.EXP; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
	  echo "OK"; \
	else \
	  echo "Error:"; \
	  diff res.out res.EXP;\
	fi
	rm -f res.out
