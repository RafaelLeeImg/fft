#!/usr/bin/python3
import math


def FFT(P):
    # coeff rep
    # P-[P(0),P(1)...P(n-1)]: P0*x^0+P1*x^1+P2*x^2+...+P(n-1)*x^(n-1)
    n = len(P)  # n is a power of 2
    if n == 1:
        return P
    k=1
    while (1<<k)<n:
        k+=1
    if (1<<k) > n:
        raise ValueError('only support array with length of 2**n')
    omega = math.e**(2*math.pi*complex(0,1)/n)
    Pe, Po = P[::2], P[1::2]
    Ye, Yo = FFT(Pe), FFT(Po)
    Y = [0]*n
    square_2 = 2**0.5
    for j in range(n//2):
        Y[j] = (Ye[j]+Yo[j]*(omega**j))/square_2
        Y[j+n//2] = (Ye[j]-Yo[j]*(omega**j))/square_2
    return Y


def IFFT(P):
    # value rep
    # P-[P(omega^0),P(omega^1)...P(omega)^(n-1)]
    n = len(P)  # n is a power of 2
    if n == 1:
        return P
    k=1
    while (1<<k)<n:
        k+=1
    if (1<<k) > n:
        raise ValueError('only support array with length of 2**n')
    omega = math.e**(-2*math.pi*complex(0,1)/n)/n
    Pe, Po = P[::2], P[1::2]
    Ye, Yo = FFT(Pe), FFT(Po)
    Y = [0]*n
    square_2 = 2**0.5
    for j in range(n//2):
        Y[j] = (Ye[j]+Yo[j]*(omega**j))/square_2
        Y[j+n//2] = (Ye[j]-Yo[j]*(omega**j))/square_2
    return Y


def main():
    # p = [2, 0, -2, 0]
    # p = [0, 2, 0, -2]
    # p = [1, 0, -1, 0]
    # p = [1, -1, 1, -1]
    # p = [2, -2, 2, -2]
    # p = [2, -2, 2, -2, 2, -2, 2, -2]
    p = [1, -1, 1, -1, 1, -1, 1, -1]
    # p = [1, -1, 1, -1, 1, -1, 1]
    # p = [1, -1]
    a = FFT(p)
    b = IFFT(a)

    print("function =", p)
    print("fft =", [round(i.real, 5) for i in a])
    print("ifft =", [round(i.real, 5) for i in b])


main()
# a = {2, -2, 2, -2, 2, -2, 2, -2};
# a = {2, -2, 2, -2};
