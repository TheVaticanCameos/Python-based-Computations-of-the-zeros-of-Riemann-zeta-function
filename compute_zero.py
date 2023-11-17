import numpy as np
from time import perf_counter
from tqdm import tqdm
from check_RH import compute_Zeta_AS, compute_Zeta_RS, FILE_NAME, CHANGE_METHOD


ACCURACY = 1e-5


def compute_zero(low, high, method, accuracy=ACCURACY):
    """
    This module computes zero of zeta function in the interval (low, high), with asserted accuracy.
    :param low: start of the interval
    :param high: end of the interval
    :param method: method ot be used to compute zero
    :param accuracy: accuracy
    :return: value of the zero
    """
    if np.sign(method(low)) * np.sign(method(high)) > 0:
        raise ValueError("Suspect there is no zero between {} and {}! Please check the interval.".format(low, high))
    else:
        err = high - low
        mid = (high + low) / 2
        low_value = method(low)
        mid_value = method(mid)
        high_value = method(high)
        while err > accuracy:
            if mid_value == 0:
                return mid
            elif np.sign(mid_value) * np.sign(low_value) < 0:
                high = mid
            elif np.sign(mid_value) * np.sign(high_value) < 0:
                low = mid

            err = high - low
            mid = (high + low) / 2
            low_value = method(low)
            mid_value = method(mid)
            high_value = method(high)

        digit = -int(np.log10(accuracy))
        return '%.{}f'.format(digit) % mid


def compute_from_file(filename, compute_num):
    """
    This module reads the previous interval which contains zeros, computes the zero, and write into file.
    :param filename: the file to read interval from
    :param compute_num: the total number of zeros to be computed
    :return: none
    """
    f_read = open(filename, 'r')
    f_write = open("./compute_zeros.txt", 'w')
    read_num = 0
    t1 = perf_counter()

    with tqdm(total=compute_num, desc='Process', leave=True, ncols=100, unit='zero', unit_scale=True) as pbar:
        while read_num < compute_num:
            line = f_read.readline()
            t = line.index('\t')
            data = line[t+2: -2]
            comma = data.index(',')
            low = float(data[0: comma])
            high = float(data[comma+1:])
            read_num += 1

            if read_num < CHANGE_METHOD:
                method = compute_Zeta_AS
            else:
                method = compute_Zeta_RS

            f_write.write("Zero No {}:\t".format(read_num))
            f_write.write(str(compute_zero(low, high, method)))
            f_write.write("\n")

            pbar.update(1)

    t2 = perf_counter()
    t_cost = t2 - t1

    f_write.write("\nWith accuracy: {}\n".format(ACCURACY))
    f_write.write("Total time cost:\t{} seconds.\n".format(t_cost))
    f_write.write("Average time cost:\t{} seconds per zero.".format(t_cost/compute_num))

    f_read.close()
    f_write.close()

    return


# compute_from_file(FILE_NAME, 500)
