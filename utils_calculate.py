from entity import Node
from math import sin


def calculate_cn(node_u, node_v):
    n_u = node_u.get_neighbor()
    n_v = node_v.get_neighbor()
    return list(set(n_u).intersection(set(n_v)))


def calculate_diff(node_u, node_v):
    n_u = node_u.get_neighbor()
    n_v = node_v.get_neighbor()
    return list(set(n_u).difference(set(n_v)))


def calculate_union(node_u, node_v):
    n_u = node_u.get_neighbor()
    n_v = node_v.get_neighbor()
    return list(set(n_u).union(n_v))


def jaccard_unweight(node_u, node_v):
    return 1 - (len(calculate_cn(node_u, node_v)) + 2) / (len(calculate_union(node_u, node_v)) + 2)




