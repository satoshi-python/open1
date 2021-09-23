
from argparse import ArgumentParser
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pyemma
import glob
import multiprocessing
import os
from msmtools.analysis import mfpt

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument("-trr", type=str, help="name of trr file")
    argparser.add_argument("-pdb", type=str, help="name of pdb file")
    argparser.add_argument("-feat", type=str, help="cordinate:CA cordinate, distance: CA distance", default = "cordinate")
    argparser.add_argument("-dimention", type=str, help="pca, tica or vamp", default = "tica")
    return argparser.parse_args()

def feat(pdb, files, feat1):
    feat = pyemma.coordinates.featurizer(pdb)
    if feat1 == "cordinate":
        feat.add_all()
    elif feat1 == "distance":
        feat.add_distances_ca()
    else:
        print("chose option -feat == cordinate or distance")
    data = pyemma.coordinates.load(files, features=feat)
    return data

def Dime(num, data):
    num_ra = 0
    if num == "pca":
        num_ra = 0
        num_dim = 500
        while True:
            pca = pyemma.coordinates.pca(data, dim=num_dim)
            output = pca.get_output()
            num_ra = cordinate(pca)
            print("ratio:", num_ra)
            if num_ra > 50:
                break
            else:
                num_dim += 200
    elif num == "tica":
        num_ra = 0
        num_dim = 500
        while True:
            tica = pyemma.coordinates.tica(data, dim=num_dim, lag=1)
            output = tica.get_output()
            num_ra = cordinate(tica)
            print("ratio:", num_ra)
            if num_ra > 50:
                break
            else:
                num_dim += 200
    elif num == "vamp":
        vamp = pyemma.coordinates.vamp(data, dim=2, lag=1)
        output = vamp.get_output()
        name = "vamp"
    else:
        print("chose option -dimention == pca, tica or vamp")
    return output

def cordinate(cor):
    num = 0
    text = cor.eigenvalues
    f = open(str(cor.dim) + "_eigenvalues.txt","w")
    for i in range(int(cor.dim)):
        num += float(text[i] / sum(text)) * 100
        a = str((text[i] / sum(text)) * 100) + " "
        f.write(a)
    f.close()
    """
    text = cor.eigenvectors
    f = open("eigenvectors.txt","w")
    for i in text:
        for j in i:
            a = str(i) + " "
            f.write(a)
        f.write("\n")
    f.close()
    """
    return num

def cluster_move_emsm(dime, name2):
    SCORE = []
    SCORE_cv = []
    dime1 = np.concatenate(dime)
    a2 = [int(i) for i in range(10,26)]
    try:
        for clusters in a2:
            cluster = pyemma.coordinates.cluster_kmeans(dime1, k=clusters, max_iter=100, stride=10, fixed_seed=1)
            dtrajs_concatenated = np.concatenate(cluster.dtrajs)
            LAG = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
            for ll in LAG:
                msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=int(ll), dt_traj='1 ps')
                MSM_vac(msm, "estimate_markov_model_cluster_number"+str(clusters)+"_LAG" + str(ll))
                SCORE.append([clusters, ll, msm.score(cluster.dtrajs)])
                SCORE_cv.append([clusters, ll, msm.score_cv(cluster.dtrajs)])
                fig = pyemma.plots.plot_markov_model(msm)
                fig[0].savefig("estimate_markov_model_move_cluster_number"+str(clusters)+"_LAG" + str(ll) + ".pdf", bbox_inches="tight")
                plt.clf()
                plt.close()
                a = msm.get_model_params()
                msm_length = len(a["P"])
                f = open("estimate_markov_model_move_cluster_number"+str(clusters)+"_LAG" + str(ll) + ".txt","w")
                a1 = "0"
                for i in range(1,msm_length):
                    a1 = a1 + " " + str(i)
                a1 = a1 + "\n"
                f.write(a1)
                for i in range(msm_length):
                    a1= str(i)
                    for j in range(msm_length):
                        a1 =a1 + " " + str(a["P"][i][j])
                    f.write(a1)
                    f.write("\n")
                f.close()
                msm.pcca(len(msm.transition_matrix))
                A = msm.metastable_sets[1]
                B = msm.metastable_sets[3]
                flux = pyemma.msm.tpt(msm, A, B)
                highest_membership = msm.metastable_distributions.argmax(1)
                msm_samples = msm.sample_by_state(50)
                coarse_state_centers = []
                for i in highest_membership:
                    coarse_state_centers.append(list(cluster.clustercenters[i])[:2])
                coarse_state_centers = np.array(coarse_state_centers)
                plot_LABEL(len(msm.transition_matrix), flux, highest_membership, coarse_state_centers, dime1, dtrajs_concatenated, name2, ll, "estimate_markov_model", msm, clusters)
                MM = "estimate_markov_model_"+ name2
                if not os.path.exists(MM):
                    os.mkdir(MM)
                label_txt(MM, msm_samples, len(msm.transition_matrix), ll)
    except ValueError:
        print("LAG:" + str(ll) + ", cluster:" + str(clusters))
        f = open("estimate_markov_model-score_" + name2 + ".txt","w")
        for j in SCORE:
            data = "cluster:" + str(j[0])
            data = data + " LAG:" + str(j[1])
            data = data + " score:" + str(j[2]) + "\n"
            f.write(data)
        f.close()        
        f = open("estimate_markov_model-score_cv_" + name2 + ".txt","w")
        for j in SCORE_cv:
            data = "cluster:" + str(j[0])
            data = data + " LAG:" + str(j[1])
            data = data + " score:" + str(j[2]) + "\n"
            f.write(data)
        f.close()
        return
    f = open("estimate_markov_model-score_" + name2 + ".txt","w")
    for j in SCORE:
        data = "cluster:" + str(j[0])
        data = data + " LAG:" + str(j[1])
        data = data + " score:" + str(j[2]) + "\n"
        f.write(data)
    f.close()
    f = open("estimate_markov_model-score_cv_" + name2 + ".txt","w")
    for j in SCORE_cv:
        data = "cluster:" + str(j[0])
        data = data + " LAG:" + str(j[1])
        data = data + " score:" + str(j[2]) + "\n"
        f.write(data)
    f.close()


def cluster_move_bmsm(dime, name2):
    SCORE = []
    SCORE_cv = []
    dime1 = np.concatenate(dime)
    a2 = [int(i) for i in range(10,26)]
    try:
        for clusters in a2:
            cluster = pyemma.coordinates.cluster_kmeans(dime1, k=clusters, max_iter=100, stride=10, fixed_seed=1)
            dtrajs_concatenated = np.concatenate(cluster.dtrajs)
            LAG = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
            for ll in LAG:
                msm = pyemma.msm.bayesian_markov_model(cluster.dtrajs, lag=int(ll), dt_traj='1 ps')
                MSM_vac(msm, "bayesian_markov_model_cluster_number"+str(clusters)+"_LAG" + str(ll))
                SCORE.append([clusters, ll, msm.score(cluster.dtrajs)])
                #SCORE_cv.append([clusters, ll, msm.score_cv(cluster.dtrajs)])
                fig = pyemma.plots.plot_markov_model(msm)
                fig[0].savefig("bayesian_markov_model_move_cluster_number"+str(clusters)+"_LAG" + str(ll) + ".pdf", bbox_inches="tight")
                plt.clf()
                plt.close()
                a = msm.get_model_params()
                msm_length = len(a["P"])
                f = open("bayesian_markov_model_move_cluster_number"+str(clusters)+"_LAG" + str(ll) + ".txt","w")
                a1 = "0"
                for i in range(1,msm_length):
                    a1 = a1 + " " + str(i)
                a1 = a1 + "\n"
                f.write(a1)
                for i in range(msm_length):
                    a1= str(i)
                    for j in range(msm_length):
                        a1 =a1 + " " + str(a["P"][i][j])
                    f.write(a1)
                    f.write("\n")
                f.close()
                msm.pcca(len(msm.transition_matrix))
                A = msm.metastable_sets[1]
                B = msm.metastable_sets[3]
                flux = pyemma.msm.tpt(msm, A, B)
                highest_membership = msm.metastable_distributions.argmax(1)
                msm_samples = msm.sample_by_state(50)
                coarse_state_centers = []
                for i in highest_membership:
                    coarse_state_centers.append(list(cluster.clustercenters[i])[:2])
                coarse_state_centers = np.array(coarse_state_centers)
                plot_LABEL(len(msm.transition_matrix), flux, highest_membership, coarse_state_centers, dime1, dtrajs_concatenated, name2, ll, "bayesian_markov_model_", msm, clusters)
                MM = "bayesian_markov_model_"+ name2
                if not os.path.exists(MM):
                    os.mkdir(MM)
                label_txt(MM, msm_samples, len(msm.transition_matrix), ll)
    except ValueError:
        print("LAG:" + str(ll) + ", cluster:" + str(clusters))
        f = open("bayesian_markov_model-score_" + name2 + ".txt","w")
        for j in SCORE:
            data = "cluster:" + str(j[0])
            data = data + " LAG:" + str(j[1])
            data = data + " score:" + str(j[2]) + "\n"
            f.write(data)
        f.close()
        f = open("bayesian_markov_model-score_cv_" + name2 + ".txt","w")
        for j in SCORE_cv:
            data = "cluster:" + str(j[0])
            data = data + " LAG:" + str(j[1])
            data = data + " score:" + str(j[2]) + "\n"
            f.write(data)
        f.close()
        return
    f = open("bayesian_markov_model-score_" + name2 + ".txt","w")
    for j in SCORE:
        data = "cluster:" + str(j[0])
        data = data + " LAG:" + str(j[1])
        data = data + " score:" + str(j[2]) + "\n"
        f.write(data)
    f.close()
    f = open("bayesian_markov_model-score_cv_" + name2 + ".txt","w")
    for j in SCORE_cv:
        data = "cluster:" + str(j[0])
        data = data + " LAG:" + str(j[1])
        data = data + " score:" + str(j[2]) + "\n"
        f.write(data)
    f.close()

def cluster_move_ehmsm(dime, name2):
    dime1 = np.concatenate(dime)
    cluster = pyemma.coordinates.cluster_kmeans(dime1, k=1000, max_iter=100, stride=1, fixed_seed=1)
    dtrajs_concatenated = np.concatenate(cluster.dtrajs)
    a2 = [int(i) for i in range(10,26)]
    try:
        for nstates in a2:
            LAG = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
            for ll in LAG:
                msm = pyemma.msm.estimate_hidden_markov_model(cluster.dtrajs, nstates, lag=int(ll), dt_traj='10 ps')
                MSM_vac(msm, "estimate_hidden_markov_model_nstates_number"+str(nstates)+"_LAG" + str(ll))
                fig = pyemma.plots.plot_markov_model(msm)
                fig[0].savefig("estimate_hidden_markov_model_move_cluster_number"+str(nstates)+"_LAG" + str(ll) + ".pdf", bbox_inches="tight")
                plt.clf()
                plt.close()
                a = msm.get_model_params()
                msm_length = len(a["P"])
                f = open("estimate_hidden_markov_model_move_cluster_number"+str(nstates)+"_LAG" + str(ll) + ".txt","w")
                a1 = "0"
                for i in range(1,msm_length):
                    a1 = a1 + " " + str(i)
                a1 = a1 + "\n"
                f.write(a1)
                for i in range(msm_length):
                    a1= str(i)
                    for j in range(msm_length):
                        a1 =a1 + " " + str(a["P"][i][j])
                    f.write(a1)
                    f.write("\n")
                f.close()
                hmm_samples = msm.sample_by_observation_probabilities(50)
                A = [0]
                B = [nstates - 1]
                flux = pyemma.msm.tpt(msm, A, B)
                highest_membership = msm.metastable_distributions.argmax(1)
                coarse_state_centers = cluster.clustercenters[msm.observable_set[highest_membership]]
                plot_LABEL(nstates, flux, highest_membership, coarse_state_centers, dime1, dtrajs_concatenated, name2, ll, "estimate_hidden_markov_model", msm)
                MM = "estimate_hidden_markov_model_"+ name2
                if not os.path.exists(MM):
                    os.mkdir(MM)
                label_txt(MM, hmm_samples, nstates, ll)
    except ValueError:
        print("LAG:" + str(ll) + ", cluster:" + str(nstates))
        return

def label_txt(MM, sample, number, ll):
    MM = MM + "/CLUSTER_" + str(number)
    if not os.path.exists(MM):
        os.mkdir(MM)
    MM = MM + "/LAG_" + str(ll)
    if not os.path.exists(MM):
        os.mkdir(MM)
    thread_list = []
    for j in range(number):
        MMM = MM + "/" + str(j + 1)
        test = Sort_SAMPLE(sample[j][:,1])
        f = open(MM+"/"+str(j)+ "cluster.txt","w")
        for v,c in test.items():
            a = str(v)
            a = a + ":"
            for aa in c:
                a =a + " " + str(aa)
            a = a + "\n"
            f.write(a)
        f.close()

def Sort_SAMPLE(sample):
    kai = {}
    a = [[int(i)] for i in range(200)]
    for i in sample:
        a[int(i / 10001)].append(int(i % 10001))
    for i in a:
        if len(i) > 1:
            WORD = str(i[0])
            for j in range(3 - len(WORD)):
                WORD = "0" + WORD
            kai[WORD] = i[1:]
    return kai


def cluster_move_bhmsm(dime, name2):
    dime1 = np.concatenate(dime)
    cluster = pyemma.coordinates.cluster_kmeans(dime1, k=1000, max_iter=100, stride=1, fixed_seed=1)
    dtrajs_concatenated = np.concatenate(cluster.dtrajs)
    a2 = [int(i) for i in range(10,26)]
    try:
        for nstates in a2:
            LAG = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
            for ll in LAG:
                msm = pyemma.msm.bayesian_hidden_markov_model(cluster.dtrajs, nstates, lag=int(ll), dt_traj='10 ps', conf=0.95)
                MSM_vac(msm, "bayesian_hidden_markov_model_nstates_number"+str(nstates)+"_LAG" + str(ll))
                fig = pyemma.plots.plot_markov_model(msm)
                fig[0].savefig("bayesian_hidden_markov_model_move_cluster_number"+str(nstates)+"_LAG" + str(ll) + ".pdf", bbox_inches="tight")
                plt.clf()
                plt.close()
                a = msm.get_model_params()
                msm_length = len(a["P"])
                f = open("bayesian_hidden_markov_model_move_cluster_number"+str(nstates)+"_LAG" + str(ll) + ".txt","w")
                a1 = "0"
                for i in range(1,msm_length):
                    a1 = a1 + " " + str(i)
                a1 = a1 + "\n"
                f.write(a1)
                for i in range(msm_length):
                    a1= str(i)
                    for j in range(msm_length):
                        a1 =a1 + " " + str(a["P"][i][j])
                    f.write(a1)
                    f.write("\n")
                f.close()
                hmm_samples = msm.sample_by_observation_probabilities(50)
                A = [0]
                B = [nstates - 1]
                flux = pyemma.msm.tpt(msm, A, B)
                highest_membership = msm.metastable_distributions.argmax(1)
                coarse_state_centers = cluster.clustercenters[msm.observable_set[highest_membership]]
                plot_LABEL(nstates, flux, highest_membership, coarse_state_centers, dime1, dtrajs_concatenated, name2, ll, "bayesian_hidden_markov_model", msm)
                MM = "bayesian_hidden_markov_model_"+ name2
                if not os.path.exists(MM):
                    os.mkdir(MM)
                label_txt(MM, hmm_samples, nstates, ll)
    except ValueError:
        print("LAG:" + str(ll) + ", cluster:" + str(nstates))
        return

def plot_LABEL(nstates, flux, highest_membership, coarse_state_centers, dime1, dtrajs_concatenated, name2, ll, name, msm, name3):
    fig, ax = plt.subplots(figsize=(10, 7))
    LABEL = [str(i) for i in range(nstates)]
    #pyemma.plots.plot_contour(dime1[:,0], dime1[:,1], flux.committor[msm.metastable_assignments[dtrajs_concatenated]], cmap='brg', ax=ax, mask=True, cbar_label=r'committor 0 $\to$ 3', alpha=0.8, zorder=-1)
    pyemma.plots.plot_contour(dime1[:,0], dime1[:,1], flux.committor[msm.metastable_assignments[msm.dtrajs_active[0]]], cmap='brg', ax=ax, mask=True, cbar_label=r'committor 0 $\to$ 3', alpha=0.8, zorder=-1)
    pyemma.plots.plot_flux(flux, coarse_state_centers, flux.stationary_distribution, ax=ax, show_committor=False, figpadding=0, show_frame=True, state_labels=LABEL, arrow_label_format='%2.e / ps')
    ax.set_xlabel(name2 + '1')
    ax.set_ylabel(name2 + '2')
    ax.set_xlim(dime1[:, 0].min(), dime1[:, 0].max())
    ax.set_ylim(dime1[:, 1].min(), dime1[:, 1].max())
    fig.savefig("move_plot_" + str(name) +"_cluster_number"+str(name3)+"_LAG" + str(ll) + ".pdf", bbox_inches="tight")
    plt.clf()
    plt.close()
    fig, ax = plt.subplots(figsize=(10, 7))
    pyemma.plots.plot_free_energy(dime1[:,0], dime1[:,1], ax=ax, legacy=False)
    pyemma.plots.plot_flux(flux, coarse_state_centers, flux.stationary_distribution, ax=ax, show_committor=False, figpadding=0, show_frame=True, state_labels=LABEL, arrow_labels = None)
    ax.set_xlabel(name2 + '1')
    ax.set_ylabel(name2 + '2')
    ax.set_xlim(dime1[:, 0].min(), dime1[:, 0].max())
    ax.set_ylim(dime1[:, 1].min(), dime1[:, 1].max())
    fig.savefig("move_plot_" + name + "_cluster_number"+str(name3)+"_LAG" + str(ll) + "_freeenergy" + ".pdf", bbox_inches="tight")
    plt.clf()
    plt.close()

def FEAT(data, name_dime):
    data = np.concatenate(data)
    fig, axes = plt.subplots(1, 2, figsize=(8, 3), sharex=True)
    pyemma.plots.plot_density(data[:,0], data[:,1], ax=axes[0], cbar=False, logscale=True)
    axes[0].set_ylabel(str(name_dime) + "1")
    axes[0].set_xlabel(str(name_dime) + "2")
    pyemma.plots.plot_free_energy(data[:,0], data[:,1], ax=axes[1], legacy=False)
    axes[1].set_ylabel(str(name_dime) + "1")
    axes[1].set_xlabel(str(name_dime) + "2")
    fig.savefig("free-energy_density.pdf", bbox_inches="tight")
    plt.clf()
    plt.close()

def MSM_vac(msm, name):
    text = msm.eigenvalues()
    f = open(name + "_eigenvalues.txt", "w")
    for t in text:
        a = str(t) + " "
        f.write(a)
    f.close()
    text = msm.eigenvectors_left()
    f = open(name + "_eigenvectors_left.txt", "w")
    for t in text:
        for i in t:
            a = str(i) + " "
            f.write(a)
        f.write("\n")
    f.close()
    text = msm.eigenvectors_right()
    f = open(name + "_eigenvectors_right.txt", "w")
    for t in text:
        for i in t:
            a = str(i) + " "
            f.write(a)
        f.write("\n")
    f.close()
    f = open(name + "_mfpt.txt", "w")
    for i in range(msm.nstates):
        text = mfpt(msm.transition_matrix, i)
        a = str(i) + " "
        f.write(a)
        for j in text:
            a = str(j) + " "
            f.write(a)
        f.write("\n")
    f.close()
    f = open(name + "_timescale.txt", "w")
    for i in msm.timescales():
        a = str(i) + " "
        f.write(a)
    f.close()


def main():
    OPTIONs = get_option()
    file1 = sorted(glob.glob(OPTIONs.trr + "*.trr"))
    print(len(file1))
    data = feat(OPTIONs.pdb, file1, OPTIONs.feat)
    data_dime =  Dime(OPTIONs.dimention, data)
    FEAT(data_dime, OPTIONs.dimention)
    cluster_move_emsm(data_dime, OPTIONs.dimention)
    #cluster_move_ehmsm(data_dime, OPTIONs.dimention)
    cluster_move_bmsm(data_dime, OPTIONs.dimention)
    #cluster_move_bhmsm(data_dime, OPTIONs.dimention)


if __name__ == '__main__':
    main()
