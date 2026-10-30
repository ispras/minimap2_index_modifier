try:
	from setuptools import setup, Extension
except ImportError:
	from distutils.core import setup
	from distutils.extension import Extension

import os, sys, platform

sys.path.append('python')

# Generated C lives under python/: prefer this checkout's quoted headers over
# an unrelated minimap.h installed alongside an external dependency.
extra_compile_args = ['-DHAVE_KALLOC', '-iquote', os.path.dirname(os.path.abspath(__file__))]
include_dirs = ["."]

if platform.machine() in ["aarch64", "arm64"]:
	include_dirs.append("sse2neon/")
	extra_compile_args.extend(['-ftree-vectorize', '-DKSW_SSE2_ONLY', '-D__SSE2__'])
else:
	extra_compile_args.append('-msse4.1') # WARNING: ancient x86_64 CPUs don't have SSE4

def readme():
	with open('python/README.rst') as f:
		return f.read()

setup(
	name = 'mappy',
	version = '2.30',
	url = 'https://github.com/lh3/minimap2',
	description = 'Minimap2 python binding',
	long_description = readme(),
	author = 'Heng Li',
	author_email = 'lh3@me.com',
	license = 'MIT',
	keywords = 'sequence-alignment',
	scripts = ['python/minimap2.py'],
	ext_modules = [Extension('mappy',
		sources = ['python/mappy.pyx', 'align.c', 'bseq.c', 'lchain.c', 'seed.c', 'format.c', 'hit.c', 'index.c', 'pe.c', 'jump.c', 'options.c',
				   'ksw2_extd2_sse.c', 'ksw2_exts2_sse.c', 'ksw2_extz2_sse.c', 'ksw2_ll_sse.c',
				   'kalloc.c', 'kthread.c', 'map.c', 'misc.c', 'sdust.c', 'sketch.c', 'esterr.c', 'splitidx.c',
				   'linked_vcf_list.c', 'mm_lift.c', 'mm_lift_idx.c', 'mm_varalt.c', 'mm_hapdp.c', 'rmin.c'],
		depends = ['minimap.h', 'bseq.h', 'kalloc.h', 'kdq.h', 'khash.h', 'kseq.h', 'ksort.h',
				   'ksw2.h', 'kthread.h', 'kvec.h', 'mmpriv.h', 'sdust.h',
				   'python/cmappy.h', 'python/cmappy.pxd', 'linked_vcf_list.h', 'mm_lift.h',
				   'mm_varalt.h', 'mm_hapdp.h', 'rmin.h'],
		extra_compile_args = extra_compile_args,
		include_dirs = include_dirs,
		libraries = ['z', 'm', 'pthread', 'hts'])],
	classifiers = [
		'Development Status :: 5 - Production/Stable',
		'License :: OSI Approved :: MIT License',
		'Operating System :: POSIX',
		'Programming Language :: C',
		'Programming Language :: Cython',
		'Programming Language :: Python :: 2.7',
		'Programming Language :: Python :: 3',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics'],
	setup_requires=["cython"])
