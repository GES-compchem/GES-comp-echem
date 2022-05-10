import os


class Config:
    def __init__(self) -> None:
        self.ncores = self.get_ncores()

    def get_ncores(self):
        try:
            return int(os.environ["OMP_NUM_THREADS"])
        except:
            return len(os.sched_getaffinity(0))
