###################################################
# Decoding black hole interior from the radiation #
# 9/19/2020 Isaac H. Kim                          #
# MIT License                                     #
###################################################
from decode import SparseDensity

n_q = 9
n_t = 2

a = SparseDensity.rand(n_q,n_t)
a.decode(0.3)
