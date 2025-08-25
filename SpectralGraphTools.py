from collections import Counter
import numpy as np
import scipy
from scipy import sparse as sp
from scipy.sparse import csgraph as csg
import networkx as nx
import matplotlib.pyplot as plt

def get_cmap(n, name='hsv'): #Function to generate colors for graphs and plots
    return plt.get_cmap(name, n)

class Graph_Spectrum:   
    def __init__(self, G, sparse_factor:float = None, complexity = False, round_eigs = 10, weighted = False, encoded = False, E = None):
        #Inputs: 
        #   - G, graph
        #   - sparse factor, determining the number of eigenvalues to calculate (default all). If sparse, k <= order(G) - 1
        #   - complex, bool to determine if eigenvectors are complex
        #   - round_eigs, int to determine number of decimal place roundings on eigval, if no rounding then round_eigs = False
        #   - weighted, bool to determine if the spectrum should consider weight
        #   - encoded, bool to determine if the spectrum should carry encoding matrices
        #   - matrix E, optional encoded values
        self.G = G
        self.k = sparse_factor
        self.complex = complexity
        self.round_eigs = round_eigs
        self.weighted = weighted
        self.encoded = encoded
        self.E = E
        
        self.spectrum: list
        self.multiplicity: dict
        
        self.main()
        
    def construct_graph_spectrum(self):        
        if self.k != None:
            csg_G = csg.reverse_cuthill_mckee(self.G)
            self.Lap_G = nx.laplacian_matrix(csg_G)
            eigval_Lap_G, eigvec_Lap_G = sp.linalg.eigsh(self.Lap_G, k = self.k)
        else:
            self.Lap_G = sp.csr_matrix(nx.laplacian_matrix(self.G)).todense()
            eigval_Lap_G, eigvec_Lap_G = scipy.linalg.eigh(self.Lap_G)
        
        eigvec_Lap_G = np.transpose(eigvec_Lap_G)
        
        if self.round_eigs != False:
            eigval_rounded = np.round(eigval_Lap_G, decimals = self.round_eigs, out = None)
            self.spectrum = list(zip(eigval_rounded, eigvec_Lap_G))
            self.spectrum = tuple(sorted(self.spectrum, key = lambda x:x[0]))
        else:
            self.spectrum = list(zip(eigval_Lap_G, eigvec_Lap_G))
            self.spectrum = tuple(sorted(self.spectrum, key = lambda x:x[0]))
        
    def calc_multiplicity(self):
        eigvals = [i[0] for i in self.spectrum]
        self.multiplicity = dict(Counter(eigvals))
        
    def main(self):
        self.construct_graph_spectrum()
        self.calc_multiplicity()
    
    def n_spectrum(self, nm_eigvals:list = None, eigval:float = None, show_graph:bool = False):
        #Input:
        #   - nm_eigvals, list of positions of eigenvalues n->m, e.g. [1, 2, 3] => plots first, second, and third entry in self.spectrum
        #   - eigval, list of eigenvalues to plot, e.g. [0, 0.5] => plots all eigenvectors with eigenvalues 0 and 0.5 from self.spectrum
        #Preffered input is nm_eigvals, so if both are given default to nm_eigvals.
        #Output:
        #   - MatPlotLib PyPlot line plot of specified eigenvectors
        
        if show_graph == True:
            fig1, (subax0, subax1) = plt.subplots(1, 2)
            fig1.suptitle('Spectra of G & Graph of G')
            nx.draw(self.G, ax = subax1)
        else:
            fig1, subax0 = plt.subplots()
            fig1.suptitle('Spectra of G')
            
        vector_dimension = self.G.order()
        subax0.set_xlim([0, self.G.order()])
        subax0.set_ylim([-1, 1])
        
        if nm_eigvals != None:
            #Plot eigvecs of eigvals n->m
            fig_colors = get_cmap(len(nm_eigvals) + 1)
            
            for i in nm_eigvals:
                a = nm_eigvals[i]
                spectrum_x = range(0, vector_dimension)
                spectrum_y = self.spectrum[a][1]
                subax0.plot(
                    spectrum_x, spectrum_y,
                    label = f'$\\lambda_{a} = {(self.spectrum[a][0])}$',
                    color = fig_colors(i),
                    lw = 7)
        else:
            if eigval == None:
                raise ValueError('No eigenvalue range or eigenvalues specified to plot!')
            
            #Check eigval is acceptable eigval.
            for i in eigval:
                if eigval[i] not in self.multiplicity.keys():
                    raise ValueError(f'Eigenvalue {eigval[i]} not an eigenvalue of Lap_G')
            
            #Generate list of indexes of eigvals to plot
            eigvals_to_plot = []
            start = 0
            for i in eigval:
                for vec_index in range(start, len(self.spectrum)):
                    if eigval[i] == self.spectrum[vec_index][0]:
                        eigvals_to_plot.append(vec_index)
                        continue
                    elif eigval[i] != self.spectrum[vec_index][0]:
                        start = vec_index + self.multiplicity[self.spectrum[vec_index][0]] - 1
                        continue
            
            fig_colors = get_cmap(len(eigvals_to_plot) + 1)
            
            #Plot eigvecs
            for i in eigvals_to_plot:
                a = eigvals_to_plot[i]         
                spectrum_x = range(0, vector_dimension)
                spectrum_y = self.spectrum[a][1]
                subax0.plot(
                    spectrum_x, spectrum_y,
                    label = f'$\\lambda_{a} = {(self.spectrum[a][0])}$',
                    color = fig_colors(i),
                    lw = 7)
                            
        if show_graph == False:
            subax0.legend(loc = 'upper right')
        
        plt.show()

class Cluster:
    #Clusters G, using the fact that the multiplicity of lambda_0 of Lap_G is the number of connected subgraphs of G
    #Inputs:
    #   - G, graph
    #   - G_Eigs, ordered array of typles of eigenvalues and eigenvectors
    #   - E, matrix of encoded data
    #Cluster methods:
    #   - Form n clusters
    #   - Form clusters of n vertexes
    #   - Form clusters of max path n 
    
    def __init__(self, G, G_Eigs, E):
        pass
    
    def main(self, cluster_method, n_cluster, m_cluster):
        pass
    
class Percolate:
    #Performs percolation on graph G, which is the removal of edges at random.
    #Inputs:
    #   - G, graph
    #   - Vertex u, starting vertex
    #   - Vertex v, neighbor of u
    #   - p, probability of percolation
    
    def __init__(self, G, u = None, v = None, p = 1.0):
        pass
    
    def main(self):
        pass
    
 