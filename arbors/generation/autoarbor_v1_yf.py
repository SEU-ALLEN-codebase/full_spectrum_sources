#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AutoArbor_M3: The method 3 of AutoArbor Algorithm
Created on Feb 13 10:22:52 2020
Last revision: March 24, 2020

@author: Hanchuan Peng
"""

print(__doc__)

import time
import os
import csv
import pandas
import numpy as np
from numpy import linalg as LA
import math
import matplotlib.pyplot as plt

from sklearn import cluster
#from sklearn.feature_extraction import image
from sklearn.cluster import spectral_clustering

#from sklearn.cluster import MiniBatchKMeans, KMeans
from sklearn.metrics.pairwise import pairwise_distances_argmin

import argparse

#some functions

def generate_tmp_swc(filename):
    #generate a temporary swc file which removes the # lines and replace the 
    #column index row
    fin = open(filename)
    fout = open(filename + "_tmp_swc.csv", "wt")
    fout.write( 'id type x y z r pid\n' )
    for line in fin:
        if not ( line.startswith('#') or (not line[0].isdigit()) ):
            ### ESWC format
            linesplit=line.split(' ')
            n=len(linesplit)
            if n>7:
                new_line=''
                for i in range(7):
                    new_line=new_line+' '+linesplit[i]
                new_line=new_line.lstrip()
                new_line=new_line+'\n'
                fout.write( new_line )
            else:
                fout.write( line )
    fin.close()
    fout.close()


def s_clustering(XS, coords, my_n_clusters):

    np.random.seed(0)
    
    # #############################################################################
    # Compute clustering 
    
    spectral = cluster.SpectralClustering(
        n_clusters=my_n_clusters, eigen_solver='arpack',random_state=0,
        affinity="precomputed")

    t0 = time.time()
    spectral.fit(XS)
    t_batch = time.time() - t0
    print(t_batch)

    my_labels = spectral.labels_
    print(my_labels)
    
    score = 0
    for k in zip(range(my_n_clusters)):
#        print(k)
        my_members = my_labels == k
        cluster_center = np.mean(coords[my_members, :], axis=0, dtype=np.float64) 
        cluster_std = np.std(coords[my_members, :], axis=0, dtype=np.float64)
#        print(cluster_std)
        score += LA.norm(cluster_std)        
 
    return [score/my_n_clusters, my_labels];    


def aarbor_adaptive_spectral_swc(filename, min_my_n_clusters, max_my_n_clusters):
    
    generate_tmp_swc(filename)
    
    tmpfile = filename + '_tmp_swc.csv'
    df = pandas.read_csv(tmpfile,
                         sep=' ')
    os.remove(f'{tmpfile}')

    X = df[['x', 'y', 'z']]
    Y = X.values

    XS = np.zeros((len(X),len(X)))

#build the index table
    xyz3 = np.zeros((len(X),3))
    for i in range(len(X)):     
        cid = df['id'][i]-1
        xyz3[cid,:] = [Y[i,0], Y[i,1], Y[i,2]]    
#        XS[i,i] = 0.999;
    
    for i in range(len(X)):     

        cp = df['pid'][i]-1 #current parent id
#        print(cp)
        
        if cp<=0:
            continue;
        
        cid = df['id'][i]-1
        dx = (xyz3[i,0]-xyz3[cp,0])
        dy = (xyz3[i,1]-xyz3[cp,1])
        dz = (xyz3[i,2]-xyz3[cp,2])
        
        XS[cid, cp] = XS[cp, cid] = math.exp(-math.sqrt(dx*dx+dy*dy+dz*dz)/100)
    #    print(X)
    
    th = 0.0001    # ln(0.0001) = -9.2
    for i in range(len(X)):     
        for j in range(len(X)):     
            if XS[i][j] < th:
                dx = (xyz3[i,0]-xyz3[j,0])
                dy = (xyz3[i,1]-xyz3[j,1])
                dz = (xyz3[i,2]-xyz3[j,2])
                XS[i][j] = math.exp(-math.sqrt(dx*dx+dy*dy+dz*dz))
                if XS[i][j]>th:
                    XS[i][j] = th;
    print(XS)
    
    # #############################################################################
    # Compute clustering 
    
    min_score = -1;
    best_n_id = 0;
    my_labels = 0;
    
    for n in range(min_my_n_clusters, max_my_n_clusters+1):
        try:
            v = s_clustering(XS, Y, n)
        except:
            return 

        cur_score = v[0]
        cur_labels = v[1]
        print([n, cur_score])
        if min_score < 0:
            min_score = cur_score;
            best_n_id = n;
            my_labels = cur_labels;            
            continue
        
        if cur_score < min_score:
            min_score = cur_score;
            best_n_id = n;
            my_labels = cur_labels;            
    
    print(my_labels)

#    my_labels = pairwise_distances_argmin(X, my_cluster_centers)

    fout = open(filename + '.autoarbor_m3.arborstat.txt' , "wt")

    
    # #############################################################################
    # Plot result
    
    fig = plt.figure(figsize=(8, 8)) #fig = plt.figure(figsize=(8, 3))
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.05, top=0.9)
    colors = ['#FF0000', '#00FF00', '#0000FF', '#FF00FF', '#FFFF00', '#00FFFF', '#F0F0F0', '#0F0F0F', '#F0FF00', '#FFF000']
    
    
    
    # KMeans
#    ax = fig.add_subplot(1, 3, 1)

    np.set_printoptions(precision=2)
    
    line = 'arbor_id arbor_node_count arbor_center_x arbor_center_y arbor_center_z\n'
    fout.write( line )
    
    ax = fig.add_subplot(1, 1, 1)
    for k, col in zip(range(best_n_id), colors):
        print(k)
        my_members = my_labels == k
        cluster_center = np.mean(Y[my_members, :], axis=0, dtype=np.float64) 
        cluster_std = np.std(Y[my_members, :], axis=0, dtype=np.float64)
        print(cluster_std)
                
        ax.plot(Y[my_members, 0], Y[my_members, 1], 'o', color=col,
                marker='o', markersize=5, label=f'arbor#{k}')
        ax.plot(cluster_center[0], cluster_center[1], '*', color=col, markerfacecolor=col,
                markeredgecolor='k', markersize=20, label=f'center of arbor#{k}')
        
        line = str(k+1) + ' ' + str(np.count_nonzero(my_members)) + ' ' + str(cluster_center)[1:-1] + '\n'
        #indeed need to calculate the arbor length in the future, 
        # not just the count of nodes
    
        fout.write( line )

    df['labels'] = my_labels #add a column  
    df.to_csv(filename + '._m3_l.eswc', index=False, sep=' ')
    df['type'] = my_labels #add a column  
    df.to_csv(filename + '._m3_lt.eswc', index=False, sep=' ')

    ax.set_title('Arbors distribution') #    ax.set_title('KMeans')
    ax.set_xticks(())
    ax.set_yticks(())
#    plt.text(-3.5, 1.8,  'train time: %.2fs\ninertia: %f' % (
#        t_batch, k_means.inertia_))
    
#    plt.show()
    plt.legend()
    plt.savefig( filename + '.autoarbor_m3.pdf')
    plt.close()
    
    fout.close()


# Main program starts here
#s=time.time()

parser = argparse.ArgumentParser()
parser.add_argument('--filename', help='SWC file name', type=str)
parser.add_argument('--L', help='Lower Bound', type=int)
parser.add_argument('--H', help='Higher Bound', type=int)
args = parser.parse_args()

aarbor_adaptive_spectral_swc(args.filename, args.L, args.H)

#e=time.time()
#t=e-s
#with open('.txt',mode='a') as f:
    #f.write(args.filename+':'+t)
    #f.write('\n')



