import random
import numpy as np




def optimizer():
    return np.random.normal(10,1)


def runBootstrap(N,alpha):
    bootstrapped_values = np.zeros(N)


    for i in range(0, N):
        bootstrapped_values[i] = optimizer()


    bootstrapped_values = np.sort(bootstrapped_values)

    low_n = int(np.ceil(N*alpha) -  1)
    high_n = int(np.ceil(N*(1-alpha)) - 1)

    low = bootstrapped_values[low_n]
    high = bootstrapped_values[high_n]

    return bootstrapped_values, low , high
