all: 
	@echo 'source("script.R")' | R --slave > nul
	@cmp -s res.out res.EXP; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
	  echo "OK"; \
	else \
	  echo "Error:"; \
	  diff res.out res.EXP;\
	fi
