import pickle
import gzip
import itertools
import random
import math
import os
import numpy as np

if __name__ == '__main__':
    # # "simulate" parameter points
    # NUM_PARAMS = 8
    # NUM_STEPS = 10
    # param_points = [
    #     [random.random() for _ in range(NUM_STEPS)]
    #     for __ in range(NUM_PARAMS)
    # ]
    #
    # # MSI configuration
    # NUM_FILES = 100 #the number of jobs we are going to run
    # NUM_CORES = 30
    # # get every "combo" of these
    # prod = itertools.product(*param_points)
    #
    #
    # batch_size = int(math.ceil(NUM_STEPS**NUM_PARAMS / NUM_FILES))
    # print('we have', batch_size / NUM_CORES, 'jobs per core')
    # chunks = itertools.batched(prod, batch_size)
    #
    #
    param_dir = 'param-files'
    os.makedirs(param_dir, exist_ok=True)
    # for i, c in enumerate(chunks):
    #     with gzip.open(f'{param_dir}/params{i}.pkl.gz', 'wb') as f:
    #         pickle.dump(c, f)
    #     if i % 20 == 0: print('done', i)


    print('example contents')
    with gzip.open(f'{param_dir}/params10.pkl.gz', 'rb') as f:
        dat = pickle.load(f)
        print(len(dat))
        print(dat[0])
        mini_dat = np.array(dat[0:10])
        split_dat = np.array_split(mini_dat, 5)
        print(mini_dat)
        print(split_dat)
        print(split_dat[0])
        # for d in dat:
        #     print(d)
