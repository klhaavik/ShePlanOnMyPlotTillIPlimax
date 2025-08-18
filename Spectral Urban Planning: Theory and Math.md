# Spectral Urban Planning: Theory and Math
## Introduction

The fundamental idea that underlies this project is that when information about urban constructs is encoded into a graph, the spectra of that graph defined through linear algebra represents certain important features of the urban construct. This document will cover that encoding process, along with how to connect concepts in urban planning theory to concepts in graph theory. One can think of the mathematical process at hand much like transforming from the time-domain to the frequency-domain when using Fourier analysis, where instead where are transforming from the spatial domain (datapoints in a geolocated database) to the spectral domain (tensors in a high-dimensional feature space). This spatial-feature transform will allow us to use the tools of linear algebra and graph theory to optimize for various outcomes and simulate various changes on the urban space being analyzed with great efficiency and with respect for the inherent interdependencies of many urban factors.

Much of this project relies on the work of urban theorists who have developed various systems, principles, morphologies, etc. which help "define the city" so to speak. While this project does not seek to engage in this part of the field, it will at times reference such frameworks to specify how mathematical analysis relates to the realities of cities and reflects the general "production" of urban space. References on urban theory which are used to develop this project include:
- A New Theory of Urban Design, by Christopher Alexander
- A Pattern Language, by Christopher Alexander
- A City is not a Tree, by Christopher Alexander
- Inventing Future Cities, by Michael Batty

Some terms used in this document are:
- Urban construct: A container for features of the urban space. E.g. building, parcel, infrastructure section, nodes in network, etc.
- Urban form: A collection of constructs which constitute a named consistuent of a city. E.g. The Chelsea Neighborhood, Green Street, residential district, etc.

## Projection of Urban Data onto Graphs
The following will be a specific presentation of graph theory in the context of urban theory. Graphs have long been used to represent urban networks, and the following will explain some of that history, as well as the unique way in which this project uses urban data to construct graphs. The graph constructed by said data will be referred to as the "projection" of the data.

A graph *G* is comprised of three sets; $V(G)$, the set of vertices; $E(G)$, the set of edges; and $\psi_G$, a function $E_i \rightarrow (V_i,V_j)$. In the most simple case, vertexes are parcels of land, and edges can represent physical adjacency of parcels. To encode more complex information, edges must touch all relevant vertexes upon which a certain urban feature depends, and a weighting tensor is used to represent the "strength" of the interaction. Many urban graphs end up being highly overlapping due to this requirement. Adjacency in "graph space" therefor is not spatial adjacnecy. A basic graph that represents this is $K_{3,3}$, known as the ["utility graph"](https://en.wikipedia.org/wiki/Planar_graph#/media/File:Biclique_K_3_3.svg) due to its common use for representing the connections of utilities like electricity and water to houses.



## Spectral Analysis

## Spatial-Feature Transform

## Projecting Graphs back onto Space
