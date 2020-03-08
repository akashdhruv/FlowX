CYTHON = cython 
CCOMP = gcc
PYINSTALLER = pyinstaller

EXE = flowx.sim
DIST = flowx.dist

PYFLAGS = --embed -3 -o 
CFLAGS = -c
EXEFLAGS = -o 

PYFILE = simulation.py
CFILE = simulation.c
OBJFILE = simulation.o

LIBDIR = /Library/Frameworks/Python.framework/Versions/3.6/include/

LIBFLAGS = -I$(LIBDIR)fno-strict-aliasing -fno-common -dynamic -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes 

$(CFILE): $(PYFILE)
	$(CYTHON) $(PYFLAGS) $(CFILE) $(PYFILE)

$(OBJFILE): $(CFILE)
	$(CCOMP) $(EXEFLAGS) $(OBJFILE) $(CFLAGS) $(CFILE) $(LIBFLAGS) `python-config --cflags`

$(EXE): $(OBJFILE)
	$(CCOMP) $(EXEFLAGS) $(EXE) $(OBJFILE) $(LIBFLAGS) `python-config --ldflags`

$(DIST) : $(PYFILE)
	$(PYINSTALLER) $(PYFILE)

clean.sim:
	rm *.c *.o $(EXE)

clean.dist:
	rm -r -f build dist simulation.spec
