import numpy as np
import matplotlib.pyplot as plt


def main():
    plt.figure(1)
    x1, y1 = np.loadtxt('sin_output.txt', usecols=(0, 1), unpack=True)
    plt.title('y = sin(x)')
    plt.xlabel('x')
    plt.ylabel('sin(x)')
    plt.plot(x1, y1)

    plt.figure(2)
    x2, y2 = np.loadtxt('cos_output.txt', usecols=(0, 1), unpack=True)
    plt.title('y = cos(x)')
    plt.xlabel('x')
    plt.ylabel('cos(x)')
    plt.plot(x2, y2)

    plt.figure(3)
    x3, y3 = np.loadtxt('atan_output.txt', usecols=(0, 1), unpack=True)
    plt.title('y = atan(x)')
    plt.xlabel('x')
    plt.ylabel('atan(x)')
    plt.plot(x3, y3)

    plt.figure(4)
    x4, y4 = np.loadtxt('exp_output.txt', usecols=(0, 1), unpack=True)
    plt.title('y = exp(x)')
    plt.xlabel('x')
    plt.ylabel('exp(x)')
    plt.plot(x4, y4)

    plt.figure(5)
    x5, y5 = np.loadtxt('log_output.txt', usecols=(0, 1), unpack=True)
    plt.title('y = log(x)')
    plt.xlabel('x')
    plt.ylabel('log(x)')
    plt.plot(x5, y5)

    plt.figure(6)
    x6, y6 = np.loadtxt('libm_sin_output.txt', usecols=(0, 1), unpack=True)
    y1_error = abs(y6 - y1)
    plt.title('sin(x) absolute error compared to libm')
    plt.xlabel('x')
    plt.ylabel('sin(x) error')
    plt.plot(x6, y1_error)
    print(f"sin maximum error : {np.nanmax(y1_error):0.60e}")

    plt.figure(7)
    x7, y7 = np.loadtxt('libm_cos_output.txt', usecols=(0, 1), unpack=True)
    y2_error = abs(y7 - y2)
    plt.title('cos(x) absolute error compared to libm')
    plt.xlabel('x')
    plt.ylabel('cos(x) error')
    plt.plot(x7, y2_error)
    print(f"cos maximum error : {np.nanmax(y2_error):0.60e}")

    plt.figure(8)
    x8, y8 = np.loadtxt('libm_atan_output.txt', usecols=(0, 1), unpack=True)
    y3_error = abs(y8 - y3)
    plt.title('atan(x) absolute error compared to libm')
    plt.xlabel('x')
    plt.ylabel('atan(x) error')
    plt.plot(x8, y3_error)
    print(f"atan maximum error : {np.nanmax(y3_error):0.60e}")

    plt.figure(9)
    x9, y9 = np.loadtxt('libm_exp_output.txt', usecols=(0, 1), unpack=True)
    y4_error = abs(y9 - y4)
    plt.title('exp(x) absolute error compared to libm')
    plt.xlabel('x')
    plt.ylabel('exp(x) error')
    plt.plot(x8, y4_error)
    print(f"exp maximum error : {np.nanmax(y4_error):0.60e}")

    plt.figure(10)
    x10, y10 = np.loadtxt('libm_log_output.txt', usecols=(0, 1), unpack=True)
    y5_error = abs(y10 - y5)
    plt.title('log(x) absolute error compared to libm')
    plt.xlabel('x')
    plt.ylabel('log(x) error')
    plt.plot(x8, y5_error)
    print(f"log maximum error : {np.nanmax(y5_error):0.60e}")

    plt.show()


if __name__ == '__main__':
    main()
