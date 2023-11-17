import numpy as np
import sympy as sym
from scipy.special import gamma, comb
import time
from tqdm import tqdm

"""
Notations: 
A complex variable s is represented as s = sigma + it.
The zeta function is the sum of n^{-s}, where n is from 1 to infinity.
The Chi function is defined as Chi(s) = 2^s * pi^{s-1} * sin(pi*s/2) * Gamma(1-s)
The Riemann-Siegel Zeta function is defined as Z(t) = Chi(1/2 + it)^{-1/2} * zeta(1/2 + it)
The error function in Riemann-Siegel method is defined as Phi_0(z) = cos(1/2*pi*z^2 + 3pi/8)/cos(pi*z),  
Phi_1(z) = 1/(12pi^2) * Phi_0'''(z)
"""

PI = np.pi
ERROR = 1e-5
CHANGE_METHOD = 200
FILE_NAME = "./zeros_records.txt"

z = sym.Symbol('z')
Phi_0 = sym.cos(0.5 * sym.pi * (z ** 2) + 3 * sym.pi / 8) / sym.cos(sym.pi * z)
Phi_1 = (1 / (12 * sym.pi * sym.pi)) * sym.diff(Phi_0, z, 3)


def compute_zeta_AS(t):

    """
    This function uses alternating series to compute zeta function.
    float:param t: imaginary part of complex variable s.
    complex:return: the value of zeta function at s = 1/2 + it.
    Note: error term is under ERROR.
    """

    mu = complex(0.5, -t)
    NUM = np.ceil((np.log2(1+2*t) + PI/2 * t * np.log2(np.e) - np.log2(ERROR) - np.log2(abs(1-2**mu))) / 3)
    # NUM is the iteration time in the computation

    def enum(k):
        """
        This function is an aid function.
        int :param k: an integer
        int :return: sum of several combinations
        """
        i = k
        sub_sum = 0
        while i <= NUM:
            sub_sum += comb(NUM, i)
            i += 1

        return sub_sum

    zeta = complex(0, 0)  # initialize the zeta value
    s = complex(0.5, t)  # complex variable s = 1/2 + it

    k = 1
    while k <= NUM:
        # sum up the first part
        zeta += (-1)**(k-1) / k**s
        k += 1

    second = 0
    while k <= 2*NUM:
        # sum up the second part
        second += (-1)**(k-1) * enum(k-NUM) / k**s
        k += 1

    second /= 2**NUM

    # get value of zeta at s
    zeta += second
    zeta /= 1-2**(1-s)

    return zeta


def compute_Zeta_AS(t):
    return (compute_zeta_AS(t) * np.exp(complex(0, compute_theta(t)))).real


def Chi(s):
    """
    This function computes the value of Chi(s).
    complex :param s: the complex variable
    complex :return: value of Chi(s)
    """
    return 2**s * PI**(s-1) * np.sin(PI*s/2) * gamma(1-s)


def compute_theta(t):
    """
    This function computes the value of theta(t) in asymptotic form.
    float :param t: imaginary part of s
    float :return: value of theta(t)
    """
    return (t / 2) * np.log(t / (2 * PI)) - (t / 2) - PI / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)


def compute_Phi(z0):
    """
    This function computes the value of Phi_0 and Phi_1.
    float :param z0: between 0 and 1
    tuple of two floats :return: (Phi_0(z0), Phi_1(z0))
    """
    return Phi_0.evalf(subs={z: z0}), Phi_1.evalf(subs={z: z0})


def compute_Zeta_RS(t):

    """
    This function uses Riemann-Siegel formula to compute Zeta function.
    float :param t: imaginary part of s(should be positive)
    float :return: value of Z(t)
    """

    tau = np.sqrt(t / (2 * PI))
    m = np.floor(tau)
    z0 = 2 * (tau - m) - 1

    Zeta = 0
    n = 1
    while n <= m:
        # sum up the main terms
        Zeta += 2 * (np.cos(compute_theta(t) - t * np.log(n))) / np.sqrt(n)
        n += 1

    # add remainders
    phi = compute_Phi(z0)
    remain = phi[0] - phi[1] / tau
    remain *= (-1) ** (m + 1)
    remain /= np.sqrt(tau)

    Zeta += remain

    return Zeta


def compute_zeta_RS(t):
    """
    This function uses Riemann-Siegel formula to compute zeta function.
    float :param t: imaginary part of s(should be positive)
    complex :return: value of zeta(1/2 + it)
    """
    mu = complex(0.5, -t)
    return np.exp(mu) * compute_Zeta_RS(t)


def zeros_numbers(T):
    """
    By analytical property, this function can compute the total number of zeros on the critical strip with 0<t<T.
    Actually, the number of zeros is the nearest integer of theta(T)/PI + 1.
    float :param T: up bound of the critical strip(should be positive)
    int :return: number of zeros in this area
    """
    estimate = compute_theta(T) / PI + 1
    if estimate - np.floor(estimate) < np.ceil(estimate) - estimate:
        return np.floor(estimate)
    else:
        return np.ceil(estimate)


def compute_Zeta(t):
    """
    This function uses different methods according to the value of T,
    in order to compute Zeta function more efficiently.
    float :param t: imaginary part of s
    float :return: value of Zeta(t)
    """
    if (t < CHANGE_METHOD) and (t > 0):
        return compute_Zeta_AS(t)
    elif t >= CHANGE_METHOD:
        return compute_Zeta_RS(t)
    else:
        raise TypeError("Argument t should be a positive real number.")


def check_RH(T, delta):
    """
    This function checks the Riemann Hypothesis with 0<t<T.
    float :param T: up bound of the discussed area(should be positive)
    float :param delta: step of sampling(should be positive)
    bool :return: True if RH is checked, else False
    """
    t1 = time.perf_counter()
    t = delta
    count_zeros = 0
    f = open(FILE_NAME, mode='w')

    with tqdm(total=np.ceil(T/delta), desc='Process', leave=True, ncols=100, unit='point', unit_scale=True) as pbar:
        while t < T:
            if np.sign(compute_Zeta(t)) * compute_Zeta(t + delta) < 0:
                count_zeros += 1
                f.write("Zero No.{}:\t({}, {})\n".format(count_zeros, np.round(t, 5), np.round(t+delta, 5)))
            t += delta
            pbar.update(1)

    print("Find {} zeros with 0<t<{}.".format(count_zeros, T))
    print("Expecting {} zeros.".format(zeros_numbers(T)))

    f.write("\n")
    f.write("Expecting {} zeros.\n".format(zeros_numbers(T)))
    f.write("Find {} zeros.\n".format(count_zeros))

    t2 = time.perf_counter()
    print("Total time cost: {} seconds.".format(t2 - t1))
    print("Average time cost: {} seconds per zero.".format(round((t2 - t1) / count_zeros, 7)))
    f.write("\nTotal time cost: {} seconds.\n".format(t2 - t1))
    f.write("Average time cost: {} seconds per zero.".format(round((t2 - t1) / count_zeros, 7)))
    f.close()

    if count_zeros == zeros_numbers(T):
        return True
    else:
        return False


# check_RH(1000, 0.1)
