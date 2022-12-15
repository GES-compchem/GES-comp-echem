import os
import logging

logger = logging.getLogger(__name__)

STRICT_MODE = True

def get_ncores():
    try:
        ncores = int(os.environ["OMP_NUM_THREADS"])
        logger.debug("Environment variable OMP_NUM_THREADS found")
    except:
        ncores = len(os.sched_getaffinity(0))
        logger.debug("Environment variable OMP_NUM_THREADS not found")

    logger.debug(f"Number of cores: {ncores}")
    return ncores

