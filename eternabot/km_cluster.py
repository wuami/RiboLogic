import argparse
import os
import settings
import numpy as np
from switch_design import read_puzzle_json
import inv_utils

def diff(secstruct1, secstruct2):
    """ returns the number of characters different between two strings """
    differences = 0
    for i in range(len(secstruct1)):
        if secstruct1[i] != secstruct2[i]:
            differences += 1
    return differences

def dist(secstructs1, secstructs2):
    """ returns the distance between two lists of strings, using sum of squares of diff"""
    distance = 0
    for i in range(len(secstructs1)):
        distance += diff(secstructs1[i], secstructs2[i])**2
    return distance

def dist_matrix(secstructs):
    """ create distance matrix """
    distance = []
    n = len(secstructs)
    for i in range(n):
        dist_row = []
        for j in range(n):
            if j < i:
                dist_row.append(distance[j][i])
            else:
                dist_row.append(dist(secstructs[i], secstructs[j]))
        distance.append(dist_row)
    return distance

def assign_clusters(secstructs, centroids):
    """ assign each secstruct list to closest centroid """
    clusters = []
    for secstruct in secstructs:
        best_distance = float("inf")
        best_index = -1
        for i, centroid in enumerate(centroids):
            distance = dist(secstruct, secstructs[centroid])
            if distance < best_distance:
                best_distance = distance
                best_index = i
        clusters.append(best_index)
    return clusters

def update_centroids(secstructs, clusters, centroids):
    """ assign centroid to secstruct with lowest max distance from other points in cluster"""
    k = len(centroids)
    new_centroids = []
    for i in range(k):
        cluster = []
        for j in range(len(clusters)):
            if clusters[j] == i:
                cluster.append(secstructs[j])
        if len(cluster) == 0:
            new_centroids.append(centroids[i])
            #new_centroids.append(np.random.randint(0,len(secstructs)))
        else:
            distances = dist_matrix(cluster)
            max_distances = np.amax(distances, axis=0)
            new_centroids.append(np.argmin(max_distances))
    return new_centroids

def silhouette(secstructs, centroids, clusters):
    distance = dist_matrix(secstructs)

    n_clusters = len(centroids)
    silh = []
    for i in range(n):
        cluster = clusters[i]
        cluster_dist = [[]]*n_clusters
        for j in range(n):
            cluster_dist[clusters[j]].append(distance[i][j])
        cluster_mean = [np.mean(x) for x in cluster_dist]
        neighbor_dist = np.amin([cluster_mean[x] for x in range(n_clusters) if x != cluster])
        silh.append((neighbor_dist-cluster_mean[cluster])/np.amax([neighbor_dist, cluster_mean[cluster]]))
    return silh

# parse arguments
p = argparse.ArgumentParser()
p.add_argument("puzzleid", help="name of puzzle filename", type=str)
args = p.parse_args()

# get sequences and secstructs
json = open(os.path.join(settings.PUZZLE_DIR, "%s.json" % args.puzzleid)).read()
puzzle = read_puzzle_json(json)
secstructs = []
with open(os.path.join(settings.PUZZLE_DIR, "%s.out" % args.puzzleid), 'r') as f:
    for line in f:
        if not line.startswith("#"):
            seq = line.split()[0]
            secstructs_current = []
            for j, target in enumerate(puzzle.targets):
                foldseq = puzzle.get_fold_sequence(seq, target)
                secstructs_current.append(inv_utils.fold(foldseq)[0])
            secstructs.append(secstructs_current)
print "read secstructs"

# kmeans cluster
k = 3
centroids = range(k)
old_centroids = range(1,k+1)
i = 0
while centroids != old_centroids:
    clusters = assign_clusters(secstructs, centroids)
    old_centroids = centroids
    centroids = update_centroids(secstructs, clusters, centroids)
    i += 1
    if i % 10 == 0:
        print "%s iterations finished" % i
        print clusters
        print centroids

print silhouette(secstructs, centroids, clusters)
