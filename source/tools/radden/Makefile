VPATH = source:../../global

radden: dqags.o precisn.o radden_modules.o
	$(COMP) $(OPT) $(OPT64I) source/radden.f90 -o $@ $^

%.o : %.f90
	$(COMP) $(OPT) $(OPT64I) -c $<

%.o : %.f
	$(COMP) $(OPT) $(OPT64I) -c $<

radden_modules.o: precisn.o
