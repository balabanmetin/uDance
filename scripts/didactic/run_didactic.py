import multiprocessing as mp
import sys

from didactic.options import *

if __name__ == "__main__":
    mp.set_start_method('fork')
    sys.setrecursionlimit(5000)

    options = options_config()

    options.func(options)

