CPP = g++ #-O3
FLAGS = -O3 -ffast-math -funroll-loops -fomit-frame-pointer -march=i686 -c #optimised for desktops
#FLAGS = -c #test
.SUFFIXES: .cpp.o
.cpp.o:
	$(CPP) $(FLAGS) $<

all:	clean NN

clean:	
	rm -f $(ALLOBJS) NN_wrap* *.o *.so

NNOBJS =	NN.o	

ALLOBJS = 	$(NNOBJS) NN

#
# To find the correct linker option
# python
# >>> import distutils.sysconfig
#>>> distutils.sysconfig.get_config_var('LINKFORSHARED')
#'-Xlinker -export-dynamic'
#


NN:	$(NNOBJS) 
	swig -python -c++ NN.i
	$(CPP) -c NN.cpp NN_wrap.cxx -I/software/Python-2.3.4/include/python2.3 -I/usr/include/python2.4
	$(CPP) -Xlinker -export-dynamic -shared NN.o NN_wrap.o  -o _NN.so 

depend:
	makedepend -I/usr/include/c++/3.2/ -- -E --  NN.cpp

# DO NOT DELETE THIS LINE -- make depend depends on it.

NN.o: NN.h /usr/include/time.h /usr/include/features.h
NN.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
NN.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs-32.h
NN.o: /usr/lib/gcc/i386-redhat-linux/4.0.2/include/stddef.h
NN.o: /usr/include/bits/time.h /usr/include/bits/types.h
NN.o: /usr/include/bits/typesizes.h /usr/include/math.h
NN.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
NN.o: /usr/include/bits/mathcalls.h
