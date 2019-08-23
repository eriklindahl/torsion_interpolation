SOURCES = geometry.c graph.c main.c interpolation.c pdbio.c rmsd.c spline.c

tinter: $(SOURCES)
	cc -O2 -o tinter $(SOURCES)

all:	tinter
	
clean:
	rm -f *o tinter
