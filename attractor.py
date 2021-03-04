from math import sin
import itertools
import networkx as nx
from entity import Node, Edge, Graph
import utils_calculate as utils
import matplotlib.pyplot as plt

path = ''
# name = 'com-amazon.ungraph.txt'
name = 'zarachy.txt'


class Attractor:
    def __init__(self, beta):
        self.beta = beta
        self.edge_num = 78
        self.node_num = 34
        self.node_list = []
        self.edge_list = []
        for i in range(self.node_num):
            node = Node(i)
            self.node_list.append(node)

    def print_graph(self):
        for i in range(len(self.edge_list)):
            print(self.edge_list[i])

    def load_data(self, path, name):
        with open(path + name, 'r') as f:
            lines = f.readlines()
            i = 0
            for line in lines:
                value = [int(s) for s in line.split()]
                a = value[0]
                b = value[1]
                edge = Edge(a - 1, b - 1, i)
                i += 1
                self.edge_list.append(edge)
                self.node_list[a - 1].add_neighbor(b - 1, i)
                self.node_list[b - 1].add_neighbor(a - 1, i)

    def find_edge(self, node_1, node_2):
        for i in range(len(self.edge_list)):
            edge = self.edge_list[i]
            if ((node_2.get_id() == edge.get_node_v()) and node_1.get_id() == edge.get_node_u()) or (
                    (node_1.get_id() == edge.get_node_v()) and node_2.get_id() == edge.get_node_u()):
                return edge.get_id()
        return None

    def node_dist(self, node_1, node_2):
        if node_1 == node_2:
            return 0
        if self.find_edge(node_1, node_2) is None:
            sum_1 = 0
            sum_2 = 0
            cn = utils.calculate_cn(node_1, node_2)
            for i in range(len(cn)):
                node_x = self.node_list[i]
                sum_1 += (self.node_dist(node_1, node_x) + self.node_dist(node_2, node_x))

            un = utils.calculate_union(node_1, node_2)
            x_y_list = list(itertools.combinations(un, 2))
            for i in range(len(x_y_list)):
                index_x = x_y_list[i][0]
                index_y = x_y_list[i][1]
                if self.find_edge(index_x, index_y) is None:
                    continue
                else:
                    sum_2 += self.node_dist(self.node_list[index_x], self.node_list[index_y])

            return 1 - sum_1 / sum_2
        else:
            return self.edge_list[self.find_edge(node_1, node_2)].get_weight()

    def DI(self, node_u, node_v):
        tep_1 = sin(1 - self.node_dist(node_u, node_v)) / node_u.get_deg()
        tep_2 = sin(1 - self.node_dist(node_u, node_v)) / node_v.get_deg()
        return 0 - (tep_1 + tep_2)

    def CI(self, node_u, node_v):
        cn = utils.calculate_cn(node_u, node_v)
        sum = 0
        for i in range(len(cn)):
            node_x = self.node_list[cn[i]]
            tep_1 = (sin(1 - self.node_dist(node_x, node_u)) * (
                    1 - self.node_dist(node_x, node_v))) / node_u.get_deg()
            tep_2 = (sin(1 - self.node_dist(node_x, node_v)) * (
                    1 - self.node_dist(node_x, node_u))) / node_v.get_deg()
            sum -= (tep_1 + tep_2)

        return sum

    def rho(self, node_x, node_v):
        # rectify d
        d = utils.jaccard_unweight(node_x, node_v)
        if 1 - d >= self.beta:
            return 1 - d
        else:
            return 1 - d - self.beta

    def EI(self, node_u, node_v):
        en_u = utils.calculate_diff(node_u, node_v)
        en_v = utils.calculate_diff(node_v, node_u)
        sum_1 = 0
        sum_2 = 0
        for i in range(len(en_u)):
            node_x = self.node_list[en_u[i]]
            sum_1 -= (sin(1 - self.node_dist(node_x, node_u)) * self.rho(node_x, node_v)) / node_u.get_deg()

        for i in range(len(en_v)):
            node_y = self.node_list[en_v[i]]
            sum_2 -= (sin(1 - self.node_dist(node_y, node_v)) * self.rho(node_y, node_u)) / node_v.get_deg()
        return sum_1 + sum_2

    def distance_init(self):
        for i in range(len(self.edge_list)):
            edge_uv = self.edge_list[i]
            node_u = self.node_list[edge_uv.get_node_u()]
            node_v = self.node_list[edge_uv.get_node_v()]
            edge_uv.set_weight(utils.jaccard_unweight(node_u, node_v))

            en_u = utils.calculate_diff(node_u, node_v)
            for j in range(len(en_u)):
                node_x = self.node_list[en_u[j]]
                edge_ux = self.edge_list[self.find_edge(node_x, node_u)]
                edge_ux.set_weight(utils.jaccard_unweight(node_u, node_x))

            en_v = utils.calculate_diff(node_v, node_u)
            for k in range(len(en_v)):
                node_y = self.node_list[en_v[k]]
                edge_vy = self.edge_list[self.find_edge(node_y, node_v)]
                edge_vy.set_weight(utils.jaccard_unweight(node_v, node_y))

    def interaction(self):
        loop = 0
        flag = True
        while flag:
            flag = False
            loop += 1
            for i in range(len(self.edge_list)):
                edge_uv = self.edge_list[i]
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

        print("loop", loop)

    def find_group(self):
        zero_list = []
        for i in range(len(self.edge_list)):
            edge = self.edge_list[i]
            if edge.get_weight() != 1:
                zero_list.append(edge.get_id())
        judge = [True for i in range(len(self.node_list))]

        group = 0
        i = 0
        while True in judge:
            if judge[i]:
                group += 1
                self.DFS(i, zero_list, judge)
                i += 1
            else:
                i += 1
        return group

    def DFS(self, node, edge_list, judge):
        judge[node] = False
        node_t = self.node_list[node]
        nei_edge = node_t.get_nei_edge()
        for i in range(len(nei_edge)):
            if nei_edge[i] in edge_list:
                edge_list.remove(nei_edge[i])
                edge = self.edge_list[nei_edge[i]]
                node_1 = edge.get_node_u()
                node_2 = edge.get_node_v()
                if judge[node_1]:
                    self.DFS(node_1, edge_list, judge)
                elif judge[node_2]:
                    self.DFS(node_2, edge_list, judge)
            else:
                continue

    def draw_network(self):
        G1 = nx.Graph()
        for i in range(len(self.node_list)):
            G1.add_node(i)
        for i in range(len(self.edge_list)):
            edge = self.edge_list[i]
            if edge.get_weight() < 1:
                G1.add_weighted_edges_from([(edge.get_node_u(), edge.get_node_v(), edge.get_weight())])

        nx.draw(G1, with_labels=True)
        plt.show()


if __name__ == '__main__':
    attractor = Attractor(0.6)
    attractor.load_data(path, name)
    attractor.distance_init()
    attractor.interaction()
    print(attractor.find_group())
