from math import sin
import networkx as nx
from entity import Node, Edge
import utils_calculate as utils
import matplotlib.pyplot as plt

path = 'data/'
name = 'football.txt'


class Attractor:
    def __init__(self, beta):
        self.beta = beta  # parameter lambda
        self.edge_num = 0  # the number of edge
        self.node_num = 0  # the number of node
        self.node_list = []  # all nodes in Graph
        self.edge_list = []  # all edges in Graph

    def init_data(self, path, name):
        with open(path + name, 'r') as f:
            lines = f.readlines()
            edges = []  # count edge and record
            for line in lines:
                val = [int(s) for s in line.split()]
                edges.append((val[0], val[1]))
                self.node_num = max(self.node_num, val[0], val[1])  # count total number of node

            # init all nodes
            for i in range(self.node_num):
                node = Node(i)
                self.node_list.append(node)

            # init Graph
            for j in range(len(edges)):
                a, b = edges[j]
                edge = Edge(a - 1, b - 1, j)
                self.edge_list.append(edge)
                self.edge_num += 1
                self.node_list[a - 1].add_neighbor(b - 1, j)
                self.node_list[b - 1].add_neighbor(a - 1, j)

    # find edge's index use node_1 and node_2
    def find_edge(self, node_1, node_2):
        if node_2.get_id() in node_1.get_neighbor():
            return node_1.get_nei_edge()[node_1.get_neighbor().index(node_2.get_id())]
        return None

    # calculate distance between node_1 and node_2
    def node_dist(self, node_1, node_2):
        # node_1 is the same as node_2
        if node_1 == node_2:
            return 0

        # the edge between node_1 and node_2 is not exist, use their cn and Jaccard
        if self.find_edge(node_1, node_2) is None:
            sum_1 = 0
            sum_2 = 0
            cn = utils.calculate_cn(node_1, node_2)
            for i in cn:
                node_x = self.node_list[i]
                sum_1 += (self.node_dist(node_1, node_x) + self.node_dist(node_2, node_x))

            # record all edges connect to node_1 or node_2
            all_edge = list(set(node_1.get_nei_edge() + node_2.get_nei_edge()))
            for j in all_edge:
                edge = self.edge_list[j]
                sum_2 += edge.get_weight()

            return 1 - sum_1 / sum_2
        else:
            # the edge between node_1 and node_2 is exist
            return self.edge_list[self.find_edge(node_1, node_2)].get_weight()

    # calculate DI
    def DI(self, node_u, node_v):
        tep_1 = sin(1 - self.node_dist(node_u, node_v)) / node_u.get_deg()
        tep_2 = sin(1 - self.node_dist(node_u, node_v)) / node_v.get_deg()
        return 0 - (tep_1 + tep_2)

    # calculate CI
    def CI(self, node_u, node_v):
        cn = utils.calculate_cn(node_u, node_v)
        sum = 0
        for i in cn:
            node_x = self.node_list[i]
            tep_1 = (sin(1 - self.node_dist(node_x, node_u)) * (
                    1 - self.node_dist(node_x, node_v))) / node_u.get_deg()
            tep_2 = (sin(1 - self.node_dist(node_x, node_v)) * (
                    1 - self.node_dist(node_x, node_u))) / node_v.get_deg()
            sum -= (tep_1 + tep_2)

        return sum

    # calculate rho, used in calculate EI
    def rho(self, node_x, node_v):
        d = self.node_dist(node_x, node_v)  # d(x,v)
        if 1 - d >= self.beta:
            return 1 - d
        else:
            return 1 - d - self.beta

    # calculate EI
    def EI(self, node_u, node_v):
        en_u = utils.calculate_diff(node_u, node_v)  # EN(u)
        en_v = utils.calculate_diff(node_v, node_u)  # EN(v)
        sum_1 = 0  # sum of x
        sum_2 = 0  # sum of y
        for i in en_u:
            node_x = self.node_list[i]
            sum_1 -= (sin(1 - self.node_dist(node_x, node_u)) * self.rho(node_x, node_v)) / node_u.get_deg()

        for i in en_v:
            node_y = self.node_list[i]
            sum_2 -= (sin(1 - self.node_dist(node_y, node_v)) * self.rho(node_y, node_u)) / node_v.get_deg()
        return sum_1 + sum_2

    # init all distance use Jaccard
    def distance_init(self):
        for edge_uv in self.edge_list:
            node_u = self.node_list[edge_uv.get_node_u()]
            node_v = self.node_list[edge_uv.get_node_v()]
            edge_uv.set_weight(utils.jaccard_unweight(node_u, node_v))

    # interaction simulate
    def interaction(self):
        loop = 0
        flag = True
        while flag:
            flag = False
            loop += 1
            for edge_uv in self.edge_list:
                if (edge_uv.get_weight() > 0) and (edge_uv.get_weight() < 1):
                    node_u = self.node_list[edge_uv.get_node_u()]
                    node_v = self.node_list[edge_uv.get_node_v()]
                    DI = self.DI(node_u, node_v)
                    CI = self.CI(node_u, node_v)
                    EI = self.EI(node_u, node_v)
                    sum_delta = DI + CI + EI
                    if sum_delta != 0:
                        weight = edge_uv.get_weight()
                        edge_uv.set_weight(weight + sum_delta)
                        if weight + sum_delta > 1:
                            edge_uv.set_weight(1)
                        if weight + sum_delta < 0:
                            edge_uv.set_weight(0)
                        flag = True
        print("loop:", loop)

    def draw_network(self):
        G1 = nx.Graph()
        for node in self.node_list:
            G1.add_node(node.get_id())
        for edge in self.edge_list:
            if edge.get_weight() < 1:
                G1.add_weighted_edges_from([(edge.get_node_u(), edge.get_node_v(), edge.get_weight())])
        print(nx.number_connected_components(G1))

    def find_group(self):
        zero_edge = []
        for edge in self.edge_list:
            if edge.get_weight() == 0:
                zero_edge.append(edge.get_id())

        visited = [0 for i in range(self.node_num)]
        group_id = 0
        while 0 in visited:
            group_id += 1
            node_id = visited.index(0)
            self.DFS(visited, node_id, group_id, zero_edge)

    def DFS(self, visited, node_id, group_id, zero_edge):
        visited[node_id] = 1
        node = self.node_list[node_id]
        node.set_group(group_id)
        nei_edge = node.get_nei_edge()
        nei = node.get_neighbor()
        i = 0
        while i < len(nei):
            if (nei_edge[i] in zero_edge) and (visited[nei[i]] == 0):
                self.DFS(visited, nei[i], group_id, zero_edge)
            else:
                i += 1

    def print_group_inf(self):
        for node in self.node_list:
            print("node_id:" + str(node.get_id()), " belongs to group_id:"+str(node.get_group()))


if __name__ == '__main__':
    attractor = Attractor(0.6)
    attractor.init_data(path, name)
    attractor.distance_init()
    attractor.interaction()
    attractor.draw_network()
    attractor.find_group()
    attractor.print_group_inf()
