class Node:
    def __init__(self, s_id):
        self.id = s_id
        self.neighbor = []
        self.nei_edge = []
        self.deg = 0
        self.group = 0

    def add_neighbor(self, n_id, e_id):
        self.neighbor.append(n_id)
        self.nei_edge.append(e_id)
        self.deg += 1

    def get_neighbor(self):
        return self.neighbor

    def get_deg(self):
        return len(self.neighbor)

    def get_id(self):
        return self.id

    def get_nei_edge(self):
        return self.nei_edge

    def get_group(self):
        return self.group

    def set_group(self, group_id):
        self.group = group_id


class Edge:
    def __init__(self, u_id, v_id, s_id):
        self.id = s_id
        self.node_u = u_id
        self.node_v = v_id
        self.weight = 1

    def get_node_u(self):
        return self.node_u

    def get_node_v(self):
        return self.node_v

    def set_weight(self, weight):
        self.weight = weight

    def get_weight(self):
        return self.weight

    def get_id(self):
        return self.id

    def __str__(self):
        return "node_u:%d, node_v:%d, weight:%f" % (self.node_u, self.node_v, self.weight)


class Graph:
    def __init__(self):
        self.node_list=[]
        self.edge_list=[]

    def get_node_list(self):
        return self.node_list

    def get_edge_list(self):
        return self.edge_list