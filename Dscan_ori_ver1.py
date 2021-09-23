
import tool_con as TC
import satoshi_cluster_ori as SC
import itertools
import multiprocessing
import MDAnalysis
from sklearn.cluster import DBSCAN

def con_thle(name, path, MED, PROB, result_list):
    u = MDAnalysis.Universe(path + "test_all.trr")
    frm = u.trajectory
    print("start:", name, "len TRR:", len(frm), "LEN PROB:", len(PROB))
    for i in range((int(name)-1) * 1000, (int(name)-1) * 1000 + 1000):
        if float(PROB[i]) >= 10 ** -10:
            cor = frm[i]
            kai = []
            kai.append(i)
            for j1, j2 in itertools.product(MED[0], MED[1]):
                kai.append(TC.con(list(cor[j1]), list(cor[j2])))
            result_list.append(kai)
    print("end:", name)

def distance(path, MED):
    manager = multiprocessing.Manager()
    result_list = manager.list()
    PROB = []
    for i in open(path + "prob.txt"):
        PROB.append(float(i))
    u = MDAnalysis.Universe(path + "test_all.trr")
    frm = u.trajectory
    NUM = [int(i) for i in range(1, int((len(frm)/1000)) + 1)]
    thread_list = []
    del frm
    for i in NUM:
        thread = multiprocessing.Process(target=con_thle,
                                         args=(i, path, MED, PROB, result_list),
                                         name=i)
        thread.start()
        thread_list.append(thread)
        if len(thread_list) == 100:
            for thread in thread_list:
                thread.join()
            thread_list = []
    print(len(result_list))
    return result_list

def read_pdb(path1):
    A = []
    B = []
    for i in open(path1):
        f = i.split()
        if f[2] == "CA":
            if f[4] == "A":
                A.append(int(f[1]) -1)
            elif f[4] == "B":
                B.append(int(f[1]) -1)
    ca = []
    ca.append(A)
    ca.append(B)
    return ca


def main():
    num_ca = read_pdb("/mnt/qnap16b/SOTU18/got/back_up/med26_aff4/pdb/HEN.pdb")
    data = distance("/mnt/qnap16b/SOTU18/got/back_up/med26_aff4/trr/",num_ca)
    for eps in range(20):
        ds = DBSCAN(eps=eps * 10 ,min_samples=100).fit(data)
        a = {}
        for one_label, frm in zip(ds.label_,data):
            if one_label != -1:
                if one_label in a:
                    a[one_label] = [(i + j)/2 for i,j in zip(a[one_label], frm)]
                else:
                    a[one_label] = list(fr,)
        f = open("aff_eps"+str(eps * 10)+"_minpts"+str(100) +".txt")
        for key, value in a.items():
            f.write("cluster:" + str(key) +"\n")
            for i in value:
                f.write(str(i) + " ")
            f.write("\n")
        f.close()


if __name__ == '__main__':
    main()
