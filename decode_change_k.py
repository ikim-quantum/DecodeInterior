###################################################
# Decoding black hole interior from the radiation #
# 9/19/2020 Isaac H. Kim                          #
# MIT License                                     #
###################################################
from decode import SparseDensity


n_q_remainder = 6
n_t_max = 6
itts = []
theta_sums = []
theta_rmss = []
wts = []

for n_t in range(n_t_max):
    n_q = n_t + n_q_remainder
    a = SparseDensity.rand(n_q,n_t)
    itt, theta_sum, theta_rms, wt_sum = a.decode_cost(0.3)
    itts.append(itt)
    theta_sums.append(theta_sum)
    theta_rmss.append(theta_rms)
    wts.append(wt_sum)

print("Iterations")
print(itts)
print("Sum of absolute value angles")
print(theta_sums)
print("RMS of angles")
print(theta_rmss)
print("Weight sum")
print(wts)
