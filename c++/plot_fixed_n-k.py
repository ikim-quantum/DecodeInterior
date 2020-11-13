import matplotlib.pyplot as plt
import numpy as np

text_n6k6 = open("n=6_k=6.txt", "r")
text_n7k6 = open("n=7_k=6.txt", "r")
text_n8k6 = open("n=8_k=6.txt", "r")
text_n9k6 = open("n=9_k=6.txt", "r")
text_n10k6 = open("n=10_k=6.txt", "r")


data_n6k6 = [float(x) for x in text_n6k6.read().split(',')]
data_n7k6 = [float(x) for x in text_n7k6.read().split(',')]
data_n8k6 = [float(x) for x in text_n8k6.read().split(',')]
data_n9k6 = [float(x) for x in text_n9k6.read().split(',')]
data_n10k6 = [float(x) for x in text_n10k6.read().split(',')]

logarithmic = True

if logarithmic:
    plt.plot(np.log(data_n6k6), label = '$n=6, n-k=6$')
    plt.plot(np.log(data_n7k6), label = '$n=7, n-k=6$')
    plt.plot(np.log(data_n8k6), label = '$n=8, n-k=6$')
    plt.plot(np.log(data_n9k6), label = '$n=9, n-k=6$')
    plt.plot(np.log(data_n10k6), label = '$n=10, n-k=6$')

    plt.xlabel("Depth")
    plt.ylabel("$\log$(Decoding complexity)")
    plt.legend(loc=2)
    plt.show()
    
else:
    plt.plot(data_n6k6, label = '$n=6, n-k=6$')
    plt.plot(data_n7k6, label = '$n=7, n-k=6$')
    plt.plot(data_n8k6, label = '$n=8, n-k=6$')
    plt.plot(data_n9k6, label = '$n=9, n-k=6$')
    plt.plot(data_n10k6, label = '$n=10, n-k=6$')

    plt.xlabel("Depth")
    plt.ylabel("$Decoding complexity")
    plt.legend(loc=2)
    plt.show()
