import os


if os.environ["OMP_NUM_THREADS"]:
    __NCORES__ = int(os.environ["OMP_NUM_THREADS"])
else:
    __NCORES__ = len(os.sched_getaffinity(0))
