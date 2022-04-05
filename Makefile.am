ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS = -I$(srcdir)/data

bin_PROGRAMS = bench_aos bench_mos bench_jastrow bench_jastrow_gpu overlap_mos

bench_aos_SOURCES = src/bench_aos.c data/Alz_small.h
bench_mos_SOURCES = src/bench_mos.c data/Alz_small.h
bench_jastrow_SOURCES = src/bench_jastrow.c data/Alz_small.h
bench_jastrow_gpu_SOURCES = src/bench_jastrow_gpu.c data/Alz_small.h
overlap_mos_SOURCES = src/overlap_mos.c data/Alz_small.h
