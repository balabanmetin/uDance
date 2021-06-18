import multiprocessing as mp

from didactic.options import *

if __name__ == "__main__":
    mp.set_start_method('fork')

    options = options_config()

    options.func(options)

