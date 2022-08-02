# -*- coding: utf-8 -*-
"""09.2 Agglomerative Hierarchical Clustering - LICA.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1UlxAHX4l4_hpvq9pi4UJoSzVEM6N8pWJ

09.2 Agglomerative Hierarchical Clustering - LICA

Author: Alexperezm | Master's End of Degree Project - 2021-2022

Objective: Based on gene expression, patients will be clustered, rising an agglomerative hierarchical clustering model, which given its tree shape, will be cut in several heights, to study if those clusters that have been developed based only on expression data, have significant differences on our query variable (Overall Survival).

Loading libraries:
"""

# Tratamiento de datos
# ==============================================================================
import numpy as np
import pandas as pd
from sklearn.datasets import make_blobs

# Gráficos
# ==============================================================================
import matplotlib.pyplot as plt
from matplotlib import style
style.use('ggplot') or plt.style.use('ggplot')

# Preprocesado y modelado
# ==============================================================================
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram
from sklearn.preprocessing import scale
from sklearn.metrics import silhouette_score
from sklearn.neighbors import kneighbors_graph
import scipy.cluster.hierarchy as sch

# Configuración warnings
# ==============================================================================
import warnings
warnings.filterwarnings('ignore')

from google.colab import drive
drive.mount('/content/drive')

"""Definition of the plot_dendogram function, self-explanatory."""

def plot_dendrogram(model, **kwargs):
    '''
    Esta función extrae la información de un modelo AgglomerativeClustering
    y representa su dendograma con la función dendogram de scipy.cluster.hierarchy
    '''
    
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)

    # Plot
    dendrogram(linkage_matrix, **kwargs)

"""Loading of the gene expression dataset:"""

df = pd.read_excel("/content/drive/MyDrive/TFM/Jessica_subread_counts.xlsx")
df.head()

"""Basic changes to adjust the dataset for the clustering procedure:"""

colnames = df["Geneid"]

colnames

PatientID=(df.columns)

del df["Geneid"]
X = df.transpose()
X.columns = colnames

X

# Escalado de datos
# ==============================================================================
X_scaled_0 = scale(X, axis=0) 
#X_scaled_1 = scale(X, axis=1)

type(X_scaled_0)

X_scaled = pd.DataFrame(X_scaled_0)
X_scaled.columns = colnames

X_scaled

X_scaled.shape

"""First of all, three different models are developed and visualized, each with different properties, being further explained in the memoir of this project; For our purpose, the "modelo_hclust_ward" will be used."""

# Modelos
# ==============================================================================
modelo_hclust_complete = AgglomerativeClustering(
                            affinity = 'euclidean',
                            linkage  = 'complete',
                            distance_threshold = 0,
                            n_clusters         = None
                        )
modelo_hclust_complete.fit(X=X_scaled)

modelo_hclust_average = AgglomerativeClustering(
                            affinity = 'euclidean',
                            linkage  = 'average',
                            distance_threshold = 0,
                            n_clusters         = None
                        )
modelo_hclust_average.fit(X=X_scaled)

modelo_hclust_ward = AgglomerativeClustering(
                            affinity = 'euclidean',
                            linkage  = 'ward',
                            distance_threshold = 0,
                            n_clusters         = None                           
                     )
modelo_hclust_ward.fit(X=X_scaled)

# Dendrogramas
# ==============================================================================
fig, axs = plt.subplots(3, 1, figsize=(8, 8))
plot_dendrogram(modelo_hclust_average, color_threshold=0, ax=axs[0])
axs[0].set_title("Distancia euclídea, Linkage average")
plot_dendrogram(modelo_hclust_complete, color_threshold=0, ax=axs[1])
axs[1].set_title("Distancia euclídea, Linkage complete")
plot_dendrogram(modelo_hclust_ward, color_threshold=0, ax=axs[2])
axs[2].set_title("Distancia euclídea, Linkage ward")
plt.tight_layout();

"""The following dendogram shows the distribution of the model that has been developed, based on the variable: "altura_corte", the user will be able to adapt the threshold to the requirements."""

fig, ax = plt.subplots(1, 1, figsize=(80, 40))
altura_corte = 500
plot_dendrogram(modelo_hclust_ward, color_threshold=altura_corte, ax=ax)
ax.set_title("Distancia euclídea, Linkage ward")
ax.axhline(y=altura_corte, c = 'black', linestyle='--', label='altura corte')
ax.legend();

"""Several results obtained from the model:"""

print(modelo_hclust_ward.n_leaves_)
print(modelo_hclust_ward.n_connected_components_)
print(modelo_hclust_ward.distances_)
print(modelo_hclust_ward.labels_)
print(modelo_hclust_ward.compute_full_tree)

modelo_hclust_ward.labels_

modelo_hclust_ward.n_clusters_

"""Now, once visualized the model and decided the optimal threshold, is time to run the clustering for the exact amount of clusters that the user will be expecting. The result will be the expected clustering, which will be studied further in R."""

modelo_hclust_ward = AgglomerativeClustering(
                            affinity = 'euclidean',
                            linkage  = 'ward',
                            #distance_threshold = 0,
                            n_clusters         = 5                           
                     )

modelo_hclust_ward.fit_predict(X_scaled)

"""Silhouette method to determine the optimal number of clusters in certain dataset, in this case it is developed for educational purposes more than its usefulness. Further explanation of the technique in the report of this project."""

# Método silhouette para identificar el número óptimo de clusters
# ==============================================================================
range_n_clusters = range(2, 20)
valores_medios_silhouette = []

for n_clusters in range_n_clusters:
    modelo = AgglomerativeClustering(
                    affinity   = 'euclidean',
                    linkage    = 'ward',
                    n_clusters = n_clusters
             )

    cluster_labels = modelo.fit_predict(X_scaled)
    silhouette_avg = silhouette_score(X_scaled, cluster_labels)
    valores_medios_silhouette.append(silhouette_avg)
    
fig, ax = plt.subplots(1, 1, figsize=(6, 3.84))
ax.plot(range_n_clusters, valores_medios_silhouette, marker='o')
ax.set_title("Evolución de media de los índices silhouette")
ax.set_xlabel('Número clusters')
ax.set_ylabel('Media índices silhouette');