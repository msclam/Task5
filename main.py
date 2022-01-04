import os
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit


def is_number(s):
    if s.count(".") == 1:  # 小数的判断
        if s[0] == "-":
            s = s[1:]
        if s[0] == ".":
            return False
        s = s.replace(".", "")
        for i in s:
            if i not in "0123456789":
                return False
        else:  # 这个else与for对应的
            return True
    elif s.count(".") == 0:  # 整数的判断
        if s[0] == "-":
            s = s[1:]
        for i in s:
            if i not in "0123456789":
                return False
        else:
            return True
    else:
        return False


def get_apophis_ra_dec(file_name):
    pos_dict = {}
    with open(file_name) as file:
        for line in file.readlines():
            line = line.strip().split()
            single_file_name = line[0][7:]
            sita = line[10]
            ra = line[20]
            dec = line[21]
            info = []
            info.append(sita)
            info.append(ra)
            info.append(dec)
            # print(single_file_name, sita, ra, dec)
            pos_dict[single_file_name] = info

    return pos_dict


def get_maxtrix_solution(file_name):
    solution_list = []
    with open(file_name) as file:
        for line in file.readlines():
            line = line.strip().split()
            if is_number(line[0]):
                solve = []
                solve.append(line[0])
                solve.append(line[1])
                # print(line[0], line[1])
                solution_list.append(solve)

    # print(solution_list)

    solution_dict = {}
    j = 0
    for i in range(357, 431, 2):
        file_name_odd = str(i) + '_2.fits'
        file_name_even = str(i + 1) + '_2.fits'
        # print(file_name_odd)
        # print(solution_list[j])
        # print(file_name_even)
        # print(solution_list[j])
        solution_dict[file_name_odd] = solution_list[j]
        solution_dict[file_name_even] = solution_list[j]
        j = j + 1

    return solution_dict


def get_coeff_a2_a1(file_name):
    coeff_dict = {}
    with open(file_name) as file:
        coeff_file_name = ""
        for line in file.readlines():
            line = line.strip().split()
            if ".txt" in line[0]:
                coeff_file_name = line[0]
                # print(coeff_file_name)
            else:
                list = []
                a2 = float(line[0])
                a1 = float(line[1])
                # print(a2, a1)
                list.append(a2)
                list.append(a1)
                coeff_dict[coeff_file_name[:3] + '_2.fits'] = list

    # print(len(coeff_dict))
    return coeff_dict


def get_new_ra_dec(pos_dict, solution_dict, coeff_dict):
    dec_dcr = {}
    for file_name in pos_dict:
        print(file_name, pos_dict[file_name])
        if file_name in solution_dict.keys() and file_name in coeff_dict.keys():
            info = pos_dict[file_name]

            sita = float(info[0])
            sita_rad = sita * np.pi / 180.0  # 度数转弧度

            bp_rp = float(solution_dict[file_name][0])

            x_tan = bp_rp * np.tan(sita_rad)

            old_y_ra_oc = float(pos_dict[file_name][1]) * 1000
            old_y_dec_oc = float(pos_dict[file_name][2]) * 1000

            a2_dec = coeff_dict[file_name][0]
            a1_dec = coeff_dict[file_name][1]

            # new_y_ra_oc = a1_ra + a2_ra * x_tan
            new_y_dec_oc = a1_dec + a2_dec * x_tan
            print("a1:{}    a2:{}   bp-rp:{}    xtan(sita):{}".format(a1_dec, a2_dec, bp_rp, x_tan))
            print("旧ra{}  旧dec{}".format(old_y_ra_oc, old_y_dec_oc))
            print("新dec: {}".format(new_y_dec_oc))
            print("deta_dec_o-c: {}".format((old_y_dec_oc - new_y_dec_oc) / 1000))
            print()

            dec_dcr[file_name] = (old_y_dec_oc - new_y_dec_oc) / 1000

    return dec_dcr

def f_1(x, a2, a1):
    return a2 * x + a1


if __name__ == "__main__":
    pos_dict = get_apophis_ra_dec('./not_j2000/veryfydata_apophis_130204.txt')
    print(pos_dict)
    print()

    solution_dict = get_maxtrix_solution('./matrix_solution/solution.txt')
    print(solution_dict)
    print()

    coeff_dict = get_coeff_a2_a1('./result_coeff/result.txt')
    print(coeff_dict)
    print()

    dec_dcr = get_new_ra_dec(pos_dict, solution_dict, coeff_dict)

    print(dec_dcr)
