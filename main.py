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
    """
    获取归算的apophis的ra dec sita
    :param file_name:
    :return:
    """
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
    """
    获取之前矩阵方程求解的 bp - rp 值
    :param file_name:
    :return:
    """
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


def get_coeff_a2_a1(file_name_dec, fille_name_ra):
    """
    通过恒星的o-c （注意ra和dec）分别导出两套a2 和 a1
    :param file_name:
    :return:
    """
    coeff_dec_dict = {}
    coeff_ra_dict = {}

    with open(file_name_dec) as file:
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
                coeff_dec_dict[coeff_file_name[:3] + '_2.fits'] = list

    with open(fille_name_ra) as file:
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
                coeff_ra_dict[coeff_file_name[:3] + '_2.fits'] = list

    # print(len(coeff_dict))
    return coeff_dec_dict, coeff_ra_dict


def get_new_ra_dec(pos_dict, solution_dict, coeff_dec_dict, coeff_ra_dict):
    """
    获得改正值new_oc
    然后返回old_oc - new_oc的结果 {文件名:改正后的结果}
    :param pos_dict:
    :param solution_dict:
    :param coeff_dict:
    :return:
    """
    dec_dcr = {}
    ra_dcr = {}
    for file_name in pos_dict:
        print(file_name, pos_dict[file_name])
        if file_name in solution_dict.keys() and file_name in coeff_dec_dict.keys() and file_name in coeff_ra_dict.keys():
            info = pos_dict[file_name]  # 旧的ra和dec的信息

            sita = float(info[0])
            sita_rad = sita * np.pi / 180.0  # 度数转弧度

            bp_rp = float(solution_dict[file_name][0])  # 两个两个方程求出的bp - rp

            x_tan = bp_rp * np.tan(sita_rad)

            old_y_ra_oc = float(pos_dict[file_name][1]) * 1000
            old_y_dec_oc = float(pos_dict[file_name][2]) * 1000

            a2_dec = coeff_dec_dict[file_name][0]
            a1_dec = coeff_dec_dict[file_name][1]

            a2_ra = coeff_ra_dict[file_name][0]
            a1_ra = coeff_ra_dict[file_name][1]

            new_y_ra_oc = a1_ra + a2_ra * x_tan
            new_y_dec_oc = a1_dec + a2_dec * x_tan  # 求出改正值

            print("a1:{}    a2:{}   bp-rp:{}    xtan(sita):{}".format(a1_dec, a2_dec, bp_rp, x_tan))
            print("旧ra{}  旧dec{}".format(old_y_ra_oc, old_y_dec_oc))
            print("新dcr的ra:{} 新dcr的dec: {}".format(new_y_ra_oc, new_y_dec_oc))
            print("deta_ra_o-c: {}".format((old_y_ra_oc - new_y_ra_oc) / 1000))
            print("deta_dec_o-c: {}".format((old_y_dec_oc - new_y_dec_oc) / 1000))
            print()

            ra_dcr[file_name] = (old_y_ra_oc - new_y_ra_oc) / 1000
            dec_dcr[file_name] = (old_y_dec_oc - new_y_dec_oc) / 1000

    return dec_dcr, ra_dcr


def get_result(ra_dcr, dec_dcr):
    with open('./reslut_ra_dec/result_ra.txt', 'w') as output_file_ra:
        for key in ra_dcr.keys():
            output_file_ra.write(key + '\t' + str(ra_dcr[key]) + '\n')

    with open('./reslut_ra_dec/result_dec.txt', 'w') as output_file_dec:
        for key in dec_dcr.keys():
            output_file_dec.write(key + '\t' + str(dec_dcr[key]) + '\n')


def f_1(x, a2, a1):
    return a2 * x + a1


if __name__ == "__main__":
    print("【归算之后的apophis的信息】")
    pos_dict = get_apophis_ra_dec('./not_j2000/veryfydata_apophis_130204.txt')
    print(pos_dict)
    print()

    print("【方程求解的bp - rp（apophis）】")
    solution_dict = get_maxtrix_solution('./matrix_solution/solution.txt')
    print(solution_dict)
    print()

    print("【恒星拟合出的a2和a1（注意Ra和Dec两个方向都需要测一次）】")
    coeff_dec_dict, coeff_ra_dict = get_coeff_a2_a1('./result_coeff/result_coeff_dec.txt',
                                                    './result_coeff/result_coeff_ra.txt')
    print(coeff_dec_dict)
    print()
    print(coeff_ra_dict)
    print()

    print("【求出DCR校正之后的ra和dec的o-c】")
    dec_dcr, ra_dcr = get_new_ra_dec(pos_dict, solution_dict, coeff_dec_dict, coeff_ra_dict)
    print()

    print("【每个文件的新ra的o-c】")
    print(ra_dcr)

    print("【每个文件的新dec的o-c】")
    print(dec_dcr)

    get_result(ra_dcr, dec_dcr)

    # res = 0
    # for key in ra_dcr.keys():
    #     res += ra_dcr[key]
    #
    # print(len(dec_dcr))
    # print(res / len(ra_dcr))
    #
    #
    # res = 0
    # with open('./not_j2000/veryfydata_apophis_130204.txt') as file:
    #     for line in file.readlines():
    #         line = line.strip().split()
    #         ra = float(line[20])
    #         dec = float(line[21])
    #         if (int(line[0][7:10]) >= 389 and int(line[0][7:10]) <= 430):
    #             res += ra
    #
    # print(res / 37)
