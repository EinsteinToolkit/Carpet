CARtest: CAR_test.F90 ../CAR_utils.F90 cctk.f cctk.h
	ifc -c -g cctk.f
	cp ../CAR_utils.F90 .
	icc -E CAR_utils.F90 > CAR_utils.f90
	ifc -c -g CAR_utils.f90
	icc -E CAR_test.F90 > CAR_test.f90
	ifc -c -g CAR_test.f90
	ifc -g -o CARtest CAR_test.o CAR_utils.o cctk.o

clean:
	rm -f *.o *.f90 *~ CARtest
