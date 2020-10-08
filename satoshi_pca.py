
from sklearn.manifold import MDS, TSNE
from sklearn.decomposition import PCA
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import pandas as pd
import seaborn as sns
import statistics
from sklearn import (manifold, datasets, decomposition, ensemble,
                     discriminant_analysis, random_projection, neighbors)

def WRITE(a1,a2,kizami,clt,path):
    print("kaki")
    f = open("{0}".format(path),"w")
    for i in range(kizami):
        for j in range(kizami):
            f.write(str(clt[j][kizami - i - 1]))
            f.write(" ")
        f.write("\n")
    f.close


def MMMMM(text):
    x = []
    x1 = []
    y1 = []
    for i in open("{0}".format(text)):
        x.append(float(i))
    for j5 in range(len(x)):
        if j5 < (len(x)/2):
            x1.append(x[j5])
        else:
            y1.append(x[j5])
    x_max = max(x1)
    x_min = min(x1)
    y_max = max(y1)
    y_min = min(y1)
    return (x_max,x_min,y_max,y_min)

def P(g,g1):
    print("Contribution of each dimension:{0}".format(np.var(g, axis=0) / np.sum(np.var(g1, axis=0))))
    print("Sum of contribution:{0}".format(np.sum(np.var(g, axis=0) / np.sum(np.var(g1, axis=0)))))
    del g,g1

def pca(vectors):
    #Determining the number of dimentions
    pca = PCA(n_components=2)
    #Data conversion for PCA
    p = pca.fit_transform(vectors)
    print("prepare PCA")
    print(pca.explained_variance_ratio_)
    P(p, vectors)
    del vectors
    return p

def P_tyuusyutu_x(text,kizami,Prob,savepath,more,mini,li="satoshi",):
    MM = MMMMM(text)
    x = []
    print(text)
    for i in open("{0}".format(text)):
        x.append(float(i))
    x1 = []
    y1 = []
    for j5 in range(len(x)):
        if j5 < (len(x)/2):
            x1.append(x[j5])
        else:
            y1.append(x[j5])
    x_max = MM[0]
    x_min = MM[1]
    y_max = MM[2]
    y_min = MM[3]
    x_jiku = (x_max - x_min)
    y_jiku = (y_max - y_min)
    x_haba = x_jiku / kizami
    y_haba = y_jiku / kizami
    prob = []
    prob_kai = []
    for i in open("{0}".format(Prob)):
        prob.append(float(i))
    for i in range(len(prob)):
        if prob[i] > 0:
            prob_kai.append(prob[i])
    clt = [[0.0 for i in range(kizami)] for j in range(kizami)]
    for j6 in range(len(x1)):
        clt[kizami - (int((y1[j6] - y_min)/y_haba) + 1)][(int((x1[j6] - x_min)/x_haba) - 1)] += float(prob_kai[j6])
    clt2 = []
    for j6 in range(kizami):
        for j7 in range(kizami):
            if clt[j6][j7] > 0.0:
                if  math.log(float(clt[j6][j7])) * (-1.986) * 0.3 < 15.0:
                    clt2.append(math.log(float(clt[j6][j7])) * (-1.986) * 0.3)
    m = float(min(clt2))
    M = float(max(clt2))
    clt2 = []
    clt5 = [[0.0 for i in range(kizami)]for j in range(kizami)]
    for j6 in range(kizami):
        for j7 in range(kizami):
            if clt[j6][j7] > 0.0:
                if math.log(float(clt[j6][j7])) * (-1.986) * 0.3 < 15.0:
                    clt5[j6][j7] += (math.log(float(clt[j6][j7])) * (-1.986) * 0.3) - m
                    clt2.append((math.log(float(clt[j6][j7])) * (-1.986) * 0.3) - m)
            else :
                clt5[j6][j7] += -1
    YY_x(clt5,kizami,y1,x1,savepath,y_min,y_haba,x_min,x_haba,more,mini)

def YY_x(clt,kizami,y1,x1,savepath,y_min,y_haba,x_min,x_haba,more,mini):    
    prob =[]
    for i in range(kizami):
        for j in range(kizami):
            prob.append(clt[i][j])
    p = []
    p = sorted(prob)
    p.reverse()
    num = [i for i in p if i <= float(more) and i > 0]
    nnum1 = [i for i in p if i <= float(more) and i >= 0 and i >= mini]
    print(len(nnum1))
    x = []
    y = []
    x2 = []
    y2 = []
    for j1 in range(len(nnum1)):
        x = []
        y = []
        for i in range(kizami):
            for j in range(kizami):
                if float(clt[i][j]) == float(nnum1[j1]):
                    x.append(int(i))
                    y.append(int(j))
        x2.append(x)
        y2.append(y)
                    
    for j in range(len(x2)):
        num = []
        num1 = nnum1[j]
        for j1 in range(len(x2[j])):
            num = []
            print(j1,"/",len(x2[j]))
            for j2 in range(len(x1)):
                if int(kizami -(int((y1[j2] - y_min)/y_haba) + 1)) == int(x2[j][j1]):
                    if int((int((x1[j2] - x_min)/x_haba) - 1)) == int(y2[j][j1]):
                        if j2 not in num:
                            num.append(j2)
                        num2 = (int((y1[j2] - y_min)/y_haba) + 1)
                        num3 = int((int((x1[j2] - x_min)/x_haba) - 1))
            kaki_x(savepath,num,num1,num2,num3)
        print("end_YY")

def kaki_x(savepath,num,num1,num2,num3):
    f = open("{0}/{1}(PC1:{3}_PC2:{2}).txt".format(savepath,str(num1),str(num2),str(num3)),"w")
    for i in range(len(num)):
        f.write(str(num[i]))
        f.write("\n")
    f.close()


def P_tyuusyutu(text,kizami,Prob,savepath,li="satoshi"):
    MM = MMMMM(text)
    x = []
    print(text)
    for i in open("{0}".format(text)):
        x.append(float(i))
    x1 = []
    y1 = []
    for j5 in range(len(x)):
        if j5 < (len(x)/2):
            x1.append(x[j5])
        else:
            y1.append(x[j5])
    x_max = MM[0]
    x_min = MM[1]
    y_max = MM[2]
    y_min = MM[3]
    x_jiku = (x_max - x_min)
    y_jiku = (y_max - y_min)
    x_haba = x_jiku / kizami
    y_haba = y_jiku / kizami
    prob = []
    prob_kai = []
    for i in open("{0}".format(Prob)):
        prob.append(float(i))
    for i in range(len(prob)):
        if prob[i] > 0:
            prob_kai.append(prob[i])
    clt = [[0.0 for i in range(kizami)] for j in range(kizami)]
    for j6 in range(len(x1)):
        clt[kizami - (int((y1[j6] - y_min)/y_haba) + 1)][(int((x1[j6] - x_min)/x_haba) - 1)] += float(prob_kai[j6])
    YY(clt,kizami,y1,x1,savepath,y_min,y_haba,x_min,x_haba)

def YY(clt,kizami,y1,x1,savepath,y_min,y_haba,x_min,x_haba):
    prob =[]
    for i in range(kizami):
        for j in range(kizami):
            prob.append(clt[i][j])
    p = []
    p = sorted(prob)
    p.reverse()
    x = []
    y = []
    x2 = []
    y2 = []
    for i in range(kizami):
        for j in range(kizami):
            if float(clt[i][j]) == float(p[0]):
                x.append(int(i))
                y.append(int(j))
                print(clt[i][j])
            if float(clt[i][j]) == float(p[1]):
                x2.append(int(i))
                y2.append(int(j))
    print("END")
    num = []
    print(float(clt[x[0]][y[0]]))
    for h2 in range(len(x)):
        for i in range(len(x1)):
            if int(kizami - (int((y1[i] - y_min)/y_haba) + 1)) == int(x[h2]):
                if int((int((x1[i] - x_min)/x_haba) - 1)) == int(y[h2]):
                    print("kizami:",kizami)
                    print(int(kizami - (int((y1[i] - y_min)/y_haba) + 1)),int((int((x1[i] - x_min)/x_haba) - 1)))
                    num.append(i)
        kaki(savepath,num)
        print("end_YY")

def kaki(savepath,num):
    f = open("{0}".format(savepath),"w")
    for i in range(len(num)):
        f.write(str(num[i]))
        f.write("\n")
    f.close()

def plot_data(data,path,a,e,ta,k1="satoshi"):
    #Preparation for drawing
    #Describe points
    #plt.scatter(data[:,0], data[:,1], c=["w" for _ in labels])
    co=0
    t = zip(data)
    x = []
    y = []
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    cp =0
    for g in t:
        if cp < int(a):
            t1  = np.array(g)
            print(t1)
            x.append(t1[0][0])
            y.append(t1[0][1])
            cp +=1
        elif cp < int(a) + int(e):
            t1  = np.array(g)
            x1.append(t1[0][0])
            y1.append(t1[0][1])
            cp +=1
        elif cp < int(ta) + int(a) + int(e):
            t1  = np.array(g)
            x2.append(t1[0][0])
            y2.append(t1[0][1])
            cp +=1
    f =  open("{0}aff.txt".format(path),"w")
    for  i in range(len(x)):
        f.write(str(x[i]))
        f.write("\n")
    for l in range(len(y)):
        f.write(str(y[l]))
        f.write("\n")
    f.close
    f = open("{0}eaf.txt".format(path),"w")
    for  i in range(len(x1)):
        f.write(str(x1[i]))
        f.write("\n")
    for l in range(len(y1)):
        f.write(str(y1[l]))
        f.write("\n")
    f.close
    f = open("{0}taf.txt".format(path),"w")
    for  i in range(len(x2)):
        f.write(str(x2[i]))
        f.write("\n")
    for l in range(len(y2)):
        f.write(str(y2[l]))
        f.write("\n")
    f.close

def PP(x,x1,x2,y,y1,y2,k1,path):
    plt.figure(figsize=(12,9))
    plt.scatter(x,y,color = "black")
    plt.savefig("{0}/eaf.png".format(path))
    plt.figure(figsize=(12,9))
    plt.scatter(x1,y1,color = "red")
    plt.savefig("{0}/aff.png".format(path))
    plt.figure(figsize=(12,9))
    plt.scatter(x2,y2,color="green")
    plt.savefig("{0}/taf.png".format(path))

def plot_datali(data,filename,labels = "satoshi",a="satoshi",k1="satoshi"):
    co=0
    t = zip(data)
    x = []
    y = []
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    cp =0
    for g in t:
        t1 = np.array(g)
        x1.append(t1[0][0])
        y1.append(t1[0][1])
    f=open("{0}".format(filename),"w")
    for i in range(len(x1)):
        f.write(str(x1[i]))
        f.write("\n")
    print(len(x1),len(y1))
    for i in range(len(y1)):
        f.write(str(y1[i]))
        f.write("\n")
    f.close()

def free_energy(path,kizami,savepath,probpath,FONT=int(14),ligand="satoshi"):
    kizami = int(kizami)
    x =[]
    x1 = []
    y1 = []
    for j3 in range(len(path)):
        print(len(path))
        x = []
        for j4 in  open("{0}".format(path[j3])):
            x.append(float(j4))
        for j5 in range(len(x)):
            if j5 < (len(x)/2):
                x1.append(x[j5])
            else:
                y1.append(x[j5])
    x_max = max(x1)
    x_min = min(x1)
    y_max = max(y1)
    y_min = min(y1)
    x_jiku = (x_max - x_min)
    y_jiku = (y_max - y_min)
    x_haba = x_jiku / kizami
    y_haba = y_jiku / kizami
    prob = []
    prob_kai = []
    for j3 in range(len(probpath)):
        for i in open("{0}".format(probpath[j3])):
            prob.append(float(i))
    for i in range(len(prob)):
        if prob[i] > 0:
            prob_kai.append(prob[i])
    x_haba = x_jiku/kizami
    y_haba = y_jiku/kizami
    clt = [[0.0 for i in range(kizami)] for j in range(kizami)]
    clt5 = [[0.0 for i in range(kizami)] for j in range(kizami)]
    x = []
    y = []
    print(len(x1))
    print(len(prob_kai))
    for i in range(kizami):
        x3 = (x_haba * i) + x_min
        x.append(x3)
        y3 = (y_haba + i) + y_min
        y.append(y3)
    for j6 in range(len(x1)):
        clt[kizami - (int((y1[j6] - y_min)/y_haba) + 1)][(int((x1[j6] - x_min)/x_haba)) -1] += float(prob_kai[j6])
    print("end")
    WRITE(ligand,"pca",kizami,clt,savepath[0])
    clt2 = []
    for j6 in range(kizami):
        for j7 in range(kizami):
            if clt[j6][j7] > 0.0:
                if  math.log(float(clt[j6][j7])) * (-1.986) * 0.3 < 15.0:
                    clt2.append(math.log(float(clt[j6][j7])) * (-1.986) * 0.3)
    m = float(min(clt2))
    M = float(max(clt2))
    clt2 = []
    for j6 in range(kizami):
        for j7 in range(kizami):
            if clt[j6][j7] > 0.0:
                if math.log(float(clt[j6][j7])) * (-1.986) * 0.3 < 15.0:
                    clt5[j6][j7] += (math.log(float(clt[j6][j7])) * (-1.986) * 0.3) - m
                    clt2.append((math.log(float(clt[j6][j7])) * (-1.986) * 0.3) - m)
            else :
                clt5[j6][j7] += (-1)
    plt.rcParams["font.size"] = 14
    df = pd.DataFrame(clt5,
                      index = [y_max - (y_haba * i) for i in range(kizami)],
                      columns = [x_min + (x_haba * j) for j in range(kizami)])
    df_mask = (df <= 0)
    cen = statistics.median(clt2) + 1
    ax = sns.heatmap(df , cmap="jet",vmin = 0 ,vmax = 10,center = cen, mask = df_mask, cbar_kws={'label': 'PMF'})
    ax.figure.axes[-1].yaxis.label.set_size(FONT)
    ax.grid()
    cax = ax.collections[0].colorbar.ax
    cax.tick_params(which='major', labelsize=FONT)
    cax.xaxis.label.set_fontsize(FONT)
    plt.xlabel("PC1", fontsize=FONT)
    plt.ylabel("PC2", fontsize=FONT)
    #plt.xticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[x_min,x_min+(x_haba * kizami/4),0,x_max-(x_haba * kizami/4),x_min])
    #plt.yticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[y_max,y_max-(y_haba * kizami/4),0,y_min+(y_haba * kizami/4),y_min])
    plt.xticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[round(x_min,2),round(x_min+(x_haba * kizami/4),2),round((x_max+x_min)/2),round(x_max-(x_haba * kizami/4),2),round(x_max,2)])
    plt.yticks([0,kizami/4,kizami/2,kizami*3/4,kizami],[round(y_max,2),round(y_max-(y_haba * kizami/4),2),round((y_max+y_min)/2),round(y_min+(y_haba * kizami/4),2),round(y_min,2)])
    plt.savefig("{0}".format(savepath[1]),bbox_inches="tight")
    plt.close()

def adjust_prob(prob_path,num_path):
    data = []
    num = []
    for i in open(num_path,"r"):
        num.append(int(i))
    num1 = 0
    num2 = 0
    for i in open(prob_path,"r"):
        if float(i) > 0:
            num1 += 1
            if num1 in num:
                data.append(float(i))
        num2 += 1
    return data 
