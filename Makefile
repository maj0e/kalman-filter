#------------------------------------------------------------------------------
# Main MakeFile for kalman project 
#------------------------------------------------------------------------------

C:
	( cd Lib ; $(MAKE) )

all: C 

library:
	( cd Lib ; $(MAKE) )

.PHONY: clean purge

clean:
	( cd Lib ; $(MAKE) clean )

purge:
	( cd Lib ; $(MAKE) purge )


