import matplotlib.pyplot as plt
import numpy as np

text_n7k4 = open("n=7_k=4.txt", "r")
text_n7k5 = open("n=7_k=5.txt", "r")
text_n7k6 = open("n=7_k=6.txt", "r")


data_n7k4 = [float(x) for x in text_n7k4.read().split(',')]
data_n7k5 = [float(x) for x in text_n7k5.read().split(',')]
data_n7k6 = [float(x) for x in text_n7k6.read().split(',')]

logarithmic = True

if logarithmic:

    plt.plot(np.log(data_n7k4), label = '$n=7, n-k=4$')
    plt.plot(np.log(data_n7k5), label = '$n=7, n-k=5$')
    plt.plot(np.log(data_n7k6), label = '$n=7, n-k=6$')

    plt.xlabel("Depth")
    plt.ylabel("$\log$(Decoding complexity)")
    plt.legend(loc=2)
    plt.show()
    
else:

    plt.plot(data_n7k4, label = '$n=7, n-k=4$')
    plt.plot(data_n7k5, label = '$n=7, n-k=5$')
    plt.plot(data_n7k6, label = '$n=7, n-k=6$')

    plt.xlabel("Depth")
    plt.ylabel("$Decoding complexity")
    plt.legend(loc=2)
    plt.show()
