import scipy.cluster
import igraph as ig
import leidenalg

def compute_genes_hierarchy_leaves(X, labels, method="complete"):
    print(X.shape)
    Z = scipy.cluster.hierarchy.linkage(X, method="complete")

    # Get the hierarchy-based ordering of genes.
    leaf_index_list = scipy.cluster.hierarchy.leaves_list(Z)
    leaf_list = labels[leaf_index_list].tolist()

    return leaf_list