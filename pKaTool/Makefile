CPP = g++ -O3
FLAGS = -c -fPIC
.SUFFIXES: .cpp.o
.cpp.o:
	$(CPP) $(FLAGS) $<

all:	clean pKa_MC 

clean:	
	rm -f $(ALLOBJS) pMC_wrap* pMC.o _pMC.so 

PKA_MCOBJS =	pMC.o

ALLOBJS = 	$(PKA_MCOBJS) 

#
# To find the correct linker option
# python
# >>> import distutils.sysconfig
#>>> distutils.sysconfig.get_config_var('LINKFORSHARED')
#'-Xlinker -export-dynamic'
#


pKa_MC:	$(PKA_MCOBJS)
	swig -python -c++ pMC.i
	$(CPP) $(FLAGS) pMC.cpp pMC_wrap.cxx -I/usr/include/python2.5 -I/usr/include/python2.6
	$(CPP) -Xlinker -export-dynamic -shared pMC.o pMC_wrap.o  -o _pMC.so 

depend:
	makedepend -I/usr/include/c++/3.2/ -- -E --  pMC.cpp 

# DO NOT DELETE THIS LINE -- make depend depends on it.

pMCo: pMC.h
