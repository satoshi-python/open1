
import numpy as np
import matplotlib.pyplot as plt
from sklearn import datasets as sd
from matplotlib.colors import ListedColormap
from sklearn.cluster import KMeans
from sklearn.cluster import MeanShift
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster import DBSCAN
from sklearn.cluster import SpectralClustering
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import estimate_bandwidth


def k_means(data,col,ran):
    feature = np.array(data)
    kmeans_model = KMeans(n_clusters=col, random_state=ran).fit(feature)
    return kmeans_model

def mean_shift(data,num_band):
    feature = np.array(data)
    #bandwidth:Note that it is less scalable than the average shift algorithm and creates a bottleneck when used.
    bandwidth = estimate_bandwidth(feature, quantile=num_band, n_samples=len(data))
    print(bandwidth)
    clustering = MeanShift(bandwidth=bandwidth).fit(feature)
    labels = clustering.labels_
    #cluster_centers = clustering.cluster_centers_
    return labels

def mean_shift_defo(data,num_band):
    feature = np.array(data)
    clustering = MeanShift(bandwidth=num_band, bin_seeding=True).fit(feature)
    labels = clustering.labels_
    return labels

def MBkmeans(data,n,bastch):
    feature = np.array(data)
    kmeans = MiniBatchKMeans(n_clusters=n, random_state=0, batch_size=bastch).fit(feature)
    return kmeans.labels_

def db_scan(data,max_point=1,mini_point=20,mini=1,maxa=30,kizami=1):
    k = []
    for K in range(mini,maxa,kizami):
        eps = K / 10
        for minPts in range(max_point,mini_point):
             dbscan = DBSCAN(eps=eps,min_samples=minPts).fit(data)
             y_dbscan = dbscan.labels_
             print("eps:",eps,",minPts:", minPts)
             # outlier数
             print(len(np.where(y_dbscan ==-1)[0]))
             # クラスタ数
             print(np.max(y_dbscan))
             # クラスタ1に含まれる点の数
             print(len(np.where(y_dbscan ==0)[0]))
             # クラスタ2に含まれる点の数
             print(len(np.where(y_dbscan ==1)[0]))
             k.append(dbscan.labels_)
    return dbscan
def db_scan_kai(data,eps,minPts):
    k = []
    eps = eps/10
    dbscan = DBSCAN(eps=abs(eps),min_samples=int(minPts)).fit(data)
    y_dbscan = dbscan.labels_
    return dbscan

def  s_cluster(data,n = 2):
    feature = np.array(data)
    clustering = SpectralClustering(n_clusters=n,assign_labels="discretize",random_state=0).fit(feature)
    return clustering.labels_

def sb_cluster(data,n= 2):
    feature = np.array(data)
    clustering = SpectralBiclustering(n_clusters=n, random_state=0).fit(feature)
    return clustering.column_labels_

def sc_cluster(data,n = 2):
    feature = np.array(data)
    clustering = SpectralCoclustering(n_clusters=n, random_state=0).fit(feature)
    return clustering.column_labels_

def AP(data):
    feature = np.array(data)
    clustering = AffinityPropagation().fit(feature)
    return clustering.labels_

def AC(data):
    feature = np.array(data)
    clustering = AgglomerativeClustering().fit(feature)
    return clustering.labels_
