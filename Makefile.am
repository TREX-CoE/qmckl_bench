ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS = -I$(srcdir)/data

bin_PROGRAMS = bench_aos bench_mos bench_jastrow


if HAVE_OMP_DEVICE
    bench_aos_SOURCES = src/bench_aos_omp.c
    bench_mos_SOURCES = src/bench_mos_omp.c
    bench_jastrow_SOURCES = src/bench_jastrow_omp.c
endif

if HAVE_ACC_DEVICE
    bench_aos_SOURCES = src/bench_aos_acc.c
    bench_mos_SOURCES = src/bench_mos_acc.c
    bench_jastrow_SOURCES = src/bench_jastrow_acc.c
endif
