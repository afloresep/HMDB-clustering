import os
import pandas as pd
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
from pyvis.network import Network
from timeit import default_timer as timer

#get input and output files
directory = os.getcwd()
tmap_frag_topology = os.path.join(directory, 'HMDB-smiles_topology.dat')
tmap_frag_coordinates = os.path.join(directory, 'HMDB-smiles_coordinates.dat')
frag_prep_tmap = os.path.join(directory, 'HMDB-smiles_prep_tmap.pkl')
centroid_file = os.path.join(directory, 'centroid.txt')

#get fragment and topology data
fragments = pd.read_pickle(frag_prep_tmap)
coordinates = list(pd.read_pickle(tmap_frag_coordinates))
coordinates = pd.DataFrame(coordinates, index=['x', 'y', 'dsmiles']).transpose()
print(coordinates.shape, flush=True)
topology = list(pd.read_pickle(tmap_frag_topology))
topology = pd.DataFrame(topology, index=['src', 'dest']).transpose()
print(topology.shape, flush=True)

##find the centroid of the TMAP tree
class Tree:
    def __init__(self, n):
        self.size = n + 1
        self.cur_size = 0
        self.tree = [[] for _ in range(self.size)]
        self.iscentroid = [False] * self.size
        self.ctree = [[] for _ in range(self.size)]

    def dfs(self, src, visited, subtree):
        visited[src] = True
        subtree[src] = 1
        self.cur_size += 1
        for adj in self.tree[src]:
            if not visited[adj] and not self.iscentroid[adj]:
                self.dfs(adj, visited, subtree)
                subtree[src] += subtree[adj]

    def findCentroid(self, src, visited, subtree):
        iscentroid = True
        visited[src] = True
        heavy_node = 0
        for adj in self.tree[src]:
            if not visited[adj] and not self.iscentroid[adj]:
                if subtree[adj] > self.cur_size//2:
                    iscentroid = False
                if heavy_node == 0 or subtree[adj] > subtree[heavy_node]:
                    heavy_node = adj
        if iscentroid and self.cur_size - subtree[src] <= self.cur_size//2:
            return src
        else:
            return self.findCentroid(heavy_node, visited, subtree)

    def findCentroidUtil(self, src):
        visited = [False] * self.size
        subtree = [0] * self.size
        self.cur_size = 0
        self.dfs(src, visited, subtree)
        for i in range(self.size):
            visited[i] = False
        centroid = self.findCentroid(src, visited, subtree)
        self.iscentroid[centroid] = True
        return centroid

    def decomposeTree(self, root):
        centroid = self.findCentroidUtil(root)
        print('centroid: ', flush=True)
        print(centroid)
        return centroid

    def addEdge(self, src, dest):
        self.tree[src].append(dest)
        self.tree[dest].append(src)

start = timer()
if __name__ == '__main__':
    tree = Tree(len(coordinates))
    for i in topology.itertuples():
         tree.addEdge(i.src, i.dest)
    centroid = tree.decomposeTree(1)

print(f'Finding centroid took: {timer() - start}sec.')

file = open(centroid_file, 'w')
file.write(str(centroid))
file.close()
