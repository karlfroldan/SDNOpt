import numpy as np 

def hit_and_miss(f, a, b, h, N):
    # h: height of the bounding box
    # N: number of trials
    x_rand = np.random.uniform(a, b, N)
    # Random y in range [0, h]
    y_rand = np.random.uniform(0, h, N)

    # (b - a) is calculated with PROD(b - a) for multidimensional

    hits = np.sum(y_rand < f(x_rand))
    estimate = (b - a) * h * (hits / N)
    return estimate

def mean_value(f, a, b, N):
    # For the domain [a ,b]
    x_rand = np.random.uniform(a, b, N)
    f_values = f(x_rand)
    f_mean = np.mean(f_values)
    estimate = f_mean * (b - a)
    return estimate

def f(x):
    return np.sin(x) ** 3

a = 0.0
b = np.pi
trials = 100000

hm = hit_and_miss(f, a, b, 1.0, trials)
mv = mean_value(f, a, b, trials)
print(hm)
print(mv)