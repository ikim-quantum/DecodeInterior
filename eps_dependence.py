###################################################
# Studying the dependence of complexity in eps    #
# 9/19/2020 Isaac H. Kim                          #
# MIT License                                     #
###################################################
from state_generator import Circuit
from decode import SparseDensity
import matplotlib.pyplot as plt
import numpy as np

n_q = 5
n_samples = 100
eps = [0.1-0.01*i for i in range(10)]

itts = []
wts = []

for i in range(len(eps)):
    print("eps={}".format(eps[i]))
    # Uncomment below if you want to study the behavior for
    # random states.
    itt_sum = 0.0
    wt_sum = 0.0
    for j in range(n_samples):
        print("eps={}, itt={}".format(i, j))
        sd = SparseDensity.rand(n_q,0)
        itt, theta1, theta2, wt = sd.decode_cost(eps[i])
        itt_sum += itt
        wt_sum += wt
        
    itts.append(itt_sum/n_samples)
    wts.append(wt_sum/n_samples)

print("itts:")
print(itts)
print("weights:")
print(wts)


xs = [np.log(1/e) for e in eps]
ys = wts
plt.plot(xs, ys)
plt.xlabel("log(1/eps)")
plt.ylabel("Number of gates")
plt.show()
