
import MDAnalysis as mda
import pyemma
import glob
import numpy as np
from argparse import ArgumentParser
import multiprocessing
import MSM_pyemma_ver3_ori as p3
import satoshi_cluster_ori as SC
import random

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument("-trr", type=str, help="name of trr file")
    argparser.add_argument("-pdb", type=str, help="name of pdb file")
    argparser.add_argument("-feat", type=str, help="cordinate:CA cordinate, distance: CA distance", default = "cordinate")
    argparser.add_argument("-dimention", type=str, help="pca, tica or vamp", default = "tica")
    return argparser.parse_args()

def get_info_trr(trr_file, info_ca, result_list, trr_num):
    print("start:", trr_file)
    u = mda.Universe(trr_file)
    frm = u.trajectory
    for num, one_frm in enumerate(frm, trr_num * len(frm)):
        trr_info = []
        trr_info.append(num)
        for ca in info_ca:
            for cordinate in range(3):
                trr_info.append(one_frm[ca][cordinate])
        result_list.append(trr_info)
    print("end:", trr_file)


def read_pdb(pdb_file):
    info_ca = []
    for i in open(pdb_file):
        f = i.split()
        if f[2] == "CA" and f[4] == "B":
            info_ca.append(int(f[1]) - 1)
    return info_ca

def read_feat():
    OPTIONs = get_option()
    files = sorted(glob.glob(OPTIONs.trr + "*.trr"))
    print("number of trr:", len(files))
    info_ca = read_pdb(OPTIONs.pdb)
    manager = multiprocessing.Manager()
    result_list = manager.list()
    thread_list = []
    if OPTIONs.feat == "cordinate":
        for num_file, one_file in enumerate(files):
            thread = multiprocessing.Process(target=get_info_trr,
                                             args=(one_file, info_ca, result_list, num_file),
                                             name=one_file)
            thread.start()
            thread_list.append(thread)
            if len(thread_list) == 100:
                for thread in thread_list:
                    thread.join()
                thread_list = []
        for thread in thread_list:
            thread.join()
        thread_list = []
    result_list = sorted(result_list, key=lambda x: x[0])
    return [f[1:] for f in result_list]

def dimensional(data):
    OPTIONs = get_option()
    return p3.Dime(OPTIONs.dimention, data)

def CLUSTER(data):
    cluster_number = 10
    return SC.k_means(data[0], cluster_number, random.randrange(100))

def main():
    OPTIONs = get_option()
    print("read_feat")
    data = read_feat()
    print("dimen")
    dim_re = dimensional(np.array(data))
    print(dim_re)
    data_clu = CLUSTER(dim_re)
    p3.cluster_move_emsm(data_clu, OPTIONs.dimention, dim_re)
    p3.cluster_move_bmsm(data_clu, OPTIONs.dimention, dim_re)

if __name__ == '__main__':
    main()
