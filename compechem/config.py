import os


try:
    __NCORES__ = int(os.environ["OMP_NUM_THREADS"])
except:
    __NCORES__ = len(os.sched_getaffinity(0))
