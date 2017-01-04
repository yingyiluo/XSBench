CC=gcc
CFLAGS=-O2 -Wall

TARGETS=librapl_reader.a etrace # sampler

all: $(TARGETS)

librapl_reader.a:  rapl_reader.o parse_cpuinfo_pkg_model.o
	rm -f $@
	ar r $@ $^

etrace:  etrace.c librapl_reader.a
	$(CC) -o $@ -O2 -Wall $< -lm -L./ -lrapl_reader 

sampler:  sampler.c librapl_reader.a
	$(CC) -o $@ -O2 -Wall $< -lm -L./ -lrapl_reader  -lpthread

clean:
	rm -f $(TARGETS)
	rm -f *.o

distclean: clean
	rm -f *~


