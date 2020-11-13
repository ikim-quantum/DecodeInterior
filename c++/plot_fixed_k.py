import matplotlib.pyplot as plt
import numpy as np

text_n4k3 = open("n=4_k=3.txt", "r")
text_n5k4 = open("n=5_k=4.txt", "r")
text_n6k5 = open("n=6_k=5.txt", "r")
text_n7k6 = open("n=7_k=6.txt", "r")


data_n4k3 = [float(x) for x in text_n4k3.read().split(',')]
data_n5k4 = [float(x) for x in text_n5k4.read().split(',')]
data_n6k5 = [float(x) for x in text_n6k5.read().split(',')]
data_n7k6 = [float(x) for x in text_n7k6.read().split(',')]

logarithmic = False

if logarithmic:

    plt.plot(np.log(data_n4k3), label = '$n=4, k=1$')
    plt.plot(np.log(data_n5k4), label = '$n=5, k=1$')
    plt.plot(np.log(data_n6k5), label = '$n=6, k=1$')
    plt.plot(np.log(data_n7k6), label = '$n=7, k=1$')

    plt.xlabel("Depth")
    plt.ylabel("$\log$(Decoding complexity)")
    plt.legend(loc=2)
    plt.show()
    
else:

    plt.plot(data_n4k3, label = '$n=4, k=1$')
    plt.plot(data_n5k4, label = '$n=5, k=1$')
    plt.plot(data_n6k5, label = '$n=6, k=1$')
    plt.plot(data_n7k6, label = '$n=7, k=1$')

    plt.xlabel("Depth")
    plt.ylabel("$Decoding complexity")
    plt.legend(loc=2)
    plt.show()
