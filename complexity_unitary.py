###################################################
# Decoding black hole interior from the radiation #
# 9/19/2020 Isaac H. Kim                          #
# MIT License                                     #
###################################################
from state_generator import Circuit
from decode import SparseDensity

n_q = 7
depth_max = 30
eps = 0.3

itts = []
wts = []

for d in range(depth_max):
    print("Depth={}".format(d))
    c = Circuit.rand_1d(n_q, d+1)
    psi = c.run()
    sd = SparseDensity()
    sd.import_from_vec(psi, 0)
    # Uncomment below if you want to study the behavior for
    # random states.
    #sd = SparseDensity.rand(n_q,0)
    itt, theta1, theta2, wt = sd.decode_cost(eps)
    itts.append(itt)
    wts.append(wt)


print("itts:")
print(itts)
print("weights:")
print(wts)
