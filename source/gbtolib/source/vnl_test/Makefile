FFLAGS = -warn all -vec-report3 -traceback -p #-O0
FC = ifort

all: test.o test_fig3.o test_fig4.o

test.o: precisn.o vnl_module.o general_quadrature.o utils.o
	$(FC) $(FFLAGS) $(OPT) test.f90 precisn.o vnl_module.o general_quadrature.o utils.o -o test

test_fig3.o: precisn.o vnl_module.o general_quadrature.o utils.o
	$(FC) $(FFLAGS) $(OPT) test_fig3.f90 precisn.o vnl_module.o general_quadrature.o utils.o -o test_fig3

test_fig4.o: precisn.o vnl_module.o general_quadrature.o utils.o
	$(FC) $(FFLAGS) $(OPT) test_fig4.f90 precisn.o vnl_module.o general_quadrature.o utils.o -o test_fig4

test: test.o test_fig3.o test_fig4.o
	./test_fig3 > out_test_fig3
	./test_fig4 > out_test_fig4
	./test > out_test

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) -c $<

clean:
	rm -f *.o *.e *~ *.mod

vnl_module.o: precisn.o general_quadrature.o const.o utils.o
utils.o: precisn.o const.o
const.o: precisn.o
general_quadrature.o: precisn.o utils.o
