import matplotlib.pyplot as plt
import numpy as np

text_n4k4 = open("n=4_k=4.txt", "r")
text_n5k5 = open("n=5_k=5.txt", "r")
text_n6k6 = open("n=6_k=6.txt", "r")


data_n4k4 = [float(x) for x in text_n4k4.read().split(',')]
data_n5k5 = [float(x) for x in text_n5k5.read().split(',')]
data_n6k6 = [float(x) for x in text_n6k6.read().split(',')]

logarithmic = False

if logarithmic:

    plt.plot(np.log(data_n4k4), label = '$n=4, n-k=4$')
    plt.plot(np.log(data_n5k5), label = '$n=5, n-k=5$')
    plt.plot(np.log(data_n6k6), label = '$n=6, n-k=6$')

    plt.xlabel("Depth")
    plt.ylabel("$\log$(Decoding complexity)")
    plt.legend(loc=2)
    plt.show()
    
else:

    plt.plot(data_n4k4, label = '$n=4, n-k=4$')
    plt.plot(data_n5k5, label = '$n=5, n-k=5$')
    plt.plot(data_n6k6, label = '$n=6, n-k=6$')

    plt.xlabel("Depth")
    plt.ylabel("$Decoding complexity")
    plt.legend(loc=2)
    plt.show()
