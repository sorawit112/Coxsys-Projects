from matplotlib import animation
from matplotlib import path
from numpy import matlib

from scipy.cluster.vq import kmeans,vq
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

import math

class infrastructureGraph:
    p = np.array([])
    con_e = np.array([])
    isWall = np.array([])
    area_edges = np.array([])
    con_a = np.array([])

    wall_graph = nx.Graph()
    area_graph = nx.Graph()
    mid_edges_graph = nx.Graph()

    wall_edges_dict = {}
    area_edges_dict = {}

    def __init__(self, selected_points, connectivity_edges, selected_wall, visual = False):
        self.p = selected_points
        self.con_e = connectivity_edges
        self.isWall = selected_wall
        
        self.VISUAL = visual

    ##############################################################################################################################
    #####                                                                                                                    #####
    #####                                             construct graph                                                        #####
    #####                                                                                                                    #####
    ##############################################################################################################################

    def initialize_graphs(self):
        self.wallGraphInit()
        self.areaGraphInit()
        self.midEdgesGraphInit()

        return None

    def wallGraphInit(self):
        name_Nodes = range(1, self.p.shape[1]+1)
        name_Edges = range(1, self.con_e[0].shape[0]+1)
        
        # create nodes
        for node in name_Nodes:
            self.wall_graph.add_node(node, Coordinates = [self.p[0][node-1], self.p[1][node-1]])

        # create edges
        for edge in name_Edges:
            w = math.hypot((self.wall_graph._node[self.con_e[0][edge-1]]['Coordinates'][0] - self.wall_graph._node[self.con_e[1][edge-1]]['Coordinates'][0]),(self.wall_graph._node[self.con_e[0][edge-1]]['Coordinates'][1] - self.wall_graph._node[self.con_e[1][edge-1]]['Coordinates'][1]))
            self.wall_edges_dict[edge] = [self.con_e[0][edge-1], self.con_e[1][edge-1]]
            self.wall_graph.add_edge(self.con_e[0][edge-1], self.con_e[1][edge-1], Number = edge, isWall = self.isWall[edge-1], midEdge = np.add([self.wall_graph._node[self.con_e[0][edge-1]]['Coordinates'][0], self.wall_graph._node[self.con_e[0][edge-1]]['Coordinates'][1]],[self.wall_graph._node[self.con_e[1][edge-1]]['Coordinates'][0], self.wall_graph._node[self.con_e[1][edge-1]]['Coordinates'][1]])/2, weight = w)

        return None
    
    def areaGraphInit(self):
        num_area = self.area_edges.shape[0]
        name_Nodes = range(1, self.p.shape[1]+1)
        vertices = []

        self.area_edges = []
        self.con_a = []

        N = 0
        explored_list = []
        all_vertices_list = []
        # automatically detect area_edges
        for node in name_Nodes:
            neighbors = [i for i in self.wall_graph.edges._adjdict[node]]
            all_vertices_list = (np.concatenate(([[1 for i in range(len(neighbors))]], [neighbors]),axis = 0).T).tolist()
            while len(all_vertices_list) != 0:
                temp_vert_list = []
                for vert in all_vertices_list:
                    before = vert[len(vert)-2]
                    neighbors = [i for i in self.wall_graph.edges._adjdict[vert[-1]] if i != before]
                
                    if len(temp_vert_list) == 0:
                        temp_vert_list = (np.concatenate((np.array([list(vert) for j in range(len(neighbors))]).T, [neighbors]), axis = 0).T).tolist()
                    else:
                        temp_vert_list = (np.concatenate((temp_vert_list, (np.concatenate((np.array([list(vert) for j in range(len(neighbors))]).T, [neighbors]), axis = 0).T).tolist()), axis = 0)).tolist()

                all_vertices_list = temp_vert_list

                for vert_list in  all_vertices_list:
                    if vert_list[0] == vert_list[-1]:
                        if len(explored_list) > 0:
                            if any([set(vert_list).issubset(set(k)) for k in explored_list]):
                                continue
                        vertices = vert_list        
                        num_ver = len(vertices)
                        self.area_edges.append([])

                        for i in range(num_ver-1):
                            self.area_edges[N].append(self.wall_graph[vertices[i]][vertices[i+1]]['Number'])
                        
                        vertices = vert_list[:-1]
                        explored_list.append(vertices) 
                        num_ver = len(vertices)
                        coordinates = [self.wall_graph._node[ver]['Coordinates'] for ver in vertices]
                        
                        centroid_i = self.centroid(coordinates)

                        ver_shift = np.subtract(coordinates, centroid_i)
                        theta_i = np.zeros((num_ver))

                        for j in range(num_ver):
                            point_a = ver_shift[j]
                            if j+1 != num_ver:
                                point_b = ver_shift[j+1]
                            else:
                                point_b = ver_shift[0]

                            theta_i[j] = np.mod(np.arctan2(point_b[1] - point_a[1], point_b[0] - point_a[0]), 2*np.pi)
                        
                        actual_theta = []
                        actual_vertices = []
                        ver_list = vertices

                        for j in range(num_ver):
                            if len(actual_theta) == 0:
                                actual_theta.append(theta_i[j])
                                actual_vertices.append(ver_list[j])
                            else:
                                if abs(theta_i[j-1] - theta_i[j]) > 0.000001:
                                    actual_theta.append(theta_i[j])
                                    actual_vertices.append(ver_list[j])

                        actual_coordinates = [self.wall_graph._node[x]['Coordinates'] for x in actual_vertices]
                        count_angle = 1
                        angle = [actual_vertices[0]]

                        for j in range(len(actual_theta)-1):
                            if abs(actual_theta[j+1] - actual_theta[j]) > 0.1:
                                count_angle += 1
                                angle.append(actual_vertices[j+1])

                        self.area_graph.add_node(N+1, centroid = self.centroid(actual_coordinates), area = PolyArea([cor[0] for cor in coordinates], [cor[1] for cor in coordinates]), vertices = vertices, actual_vertices = actual_vertices, area_edges = self.area_edges[N], angle = angle)
                    
                        N += 1                        

                deleted = []
                for i in range(len(all_vertices_list)):
                    vertices = all_vertices_list[i]
                    num_ver = len(vertices)
                    coordinates = [self.wall_graph._node[ver]['Coordinates'] for ver in vertices]
                    
                    centroid_i = self.centroid(coordinates)

                    ver_shift = np.subtract(coordinates, centroid_i)
                    theta_i = np.zeros((num_ver))

                    for j in range(num_ver):
                        point_a = ver_shift[j]
                        if j+1 != num_ver:
                            point_b = ver_shift[j+1]
                        else:
                            point_b = ver_shift[0]

                        theta_i[j] = np.mod(np.arctan2(point_b[1] - point_a[1], point_b[0] - point_a[0]), 2*np.pi)
                    
                    actual_theta = []
                    actual_vertices = []
                    ver_list = vertices

                    for j in range(num_ver):
                        if len(actual_theta) == 0:
                            actual_theta.append(theta_i[j])
                            actual_vertices.append(ver_list[j])
                        else:
                            if abs(theta_i[j-1] - theta_i[j]) > 0.000001:
                                actual_theta.append(theta_i[j])
                                actual_vertices.append(ver_list[j])

                    actual_coordinates = [self.wall_graph._node[x]['Coordinates'] for x in actual_vertices]
                    count_angle = 1

                    for j in range(len(actual_theta)-1):
                        if abs(actual_theta[j+1] - actual_theta[j]) > 0.1:
                            count_angle += 1
                            
                    if count_angle > 5:
                        deleted.append(i)
                        continue

                    for j in range(len(explored_list)):
                        angle_j = self.area_graph._node[j+1]['angle']
                        sub_list = angle_j[2] 

                        if set(explored_list[j]).issubset(set(all_vertices_list[i])):
                            deleted.append(i)
                            
                count = 0
                for i in deleted:
                    all_vertices_list.pop(i - count)
                    count += 1

        # automatically detect vertices of an area based on edge
        for i in range(num_area):
            area_edges_i = self.area_edges[i]
            wall_i = [self.wall_edges_dict[edge] for edge in area_edges_i]
            vertices_i = []

            for j in range(len(area_edges_i)):
                if j == 0:
                    vertices_i.append(list(set(wall_i[len(wall_i)-1]) & set(wall_i[j]))[0])
                else:
                    vertices_i.append(list(set(wall_i[j-1]) & set(wall_i[j]))[0])
            
            vertices.append(vertices_i)

        name_Nodes = range(1, num_area+1)
        num_edges = self.con_a.shape[1]
    
        # create nodes
        for node in name_Nodes:
            self.area_graph.add_node(node, centroid = [0, 0], area = 0, vertices = vertices[node-1], actual_vertices = [], area_edges = self.area_edges[node-1])

        # create edges
        for edge in range(num_edges):
            self.area_edges_dict[edge] = [self.con_a[0][edge], self.con_a[1][edge]]
            self.area_graph.add_edge(self.con_a[0][edge], self.con_a[1][edge])
    
        # find centroid and actual vertices 
        for i in self.area_graph._node:
            coordinates = [self.wall_graph._node[ver]['Coordinates'] for ver in self.area_graph._node[i]['vertices']]
            self.area_graph._node[i]['area'] = PolyArea([cor[0] for cor in coordinates], [cor[1] for cor in coordinates])
            num_ver = len(self.area_graph._node[i]['vertices'])
            centroid_i = self.centroid(coordinates)
            self.area_graph._node[i]['centroid'] = centroid_i

            ver_shift = np.subtract(coordinates, centroid_i)
            theta_i = np.zeros((num_ver))

            for j in range(num_ver):
                point_a = ver_shift[j]
                if j+1 != num_ver:
                    point_b = ver_shift[j+1]
                else:
                    point_b = ver_shift[0]

                theta_i[j] = np.mod(np.arctan2(point_b[1] - point_a[1], point_b[0] - point_a[0]), 2*np.pi)
            
            actual_theta = []
            actual_vertices = []
            ver_list = self.area_graph._node[i]['vertices']

            for j in range(num_ver):
                if len(actual_theta) == 0:
                    actual_theta.append(theta_i[j])
                    actual_vertices.append(ver_list[j])
                else:
                    if abs(theta_i[j-1] - theta_i[j]) > 0.000001:
                        actual_theta.append(theta_i[j])
                        actual_vertices.append(ver_list[j])

            self.area_graph._node[i]['actual_vertices'] = actual_vertices
            actual_coordinates = [self.wall_graph._node[x]['Coordinates'] for x in actual_vertices]
            self.area_graph._node[i]['centroid'] = self.centroid(actual_coordinates)

        return None

    def midEdgesGraphInit(self):
        con_mid = []

        for i in self.area_graph._node:
            edge_i = self.area_graph._node[i]['area_edges']
            isConnected = [ not self.wall_graph[self.wall_edges_dict[i][0]][self.wall_edges_dict[i][1]]['isWall'] for i in edge_i ]
            if sum(isConnected) > 1:
                temp_list = [edge_i[x] for x in range(len(edge_i)) if isConnected[x]]
                n_list = len(temp_list)
                for j in range(n_list):
                    for k in range(j+1, n_list):
                        con_mid.append([temp_list[j], temp_list[k]])

        sorted_list = [sorted(x) for x in con_mid]
        con_mid = sorted(sorted_list)
        node_mid_edge = sorted(list(set(np.array(con_mid).reshape(1,len(con_mid[0])*len(con_mid))[0])))

        # create nodes
        for i in range(len(node_mid_edge)):
            self.mid_edges_graph.add_node(i+1, edge_Number = node_mid_edge[i])

        # create edges
        for x in con_mid:
            edge_a = self.wall_edges_dict[x[0]]
            edge_b = self.wall_edges_dict[x[1]]
            point_a = self.wall_graph[edge_a[0]][edge_a[1]]['midEdge']
            point_b = self.wall_graph[edge_b[0]][edge_b[1]]['midEdge']

            mid_edge_idx = [node_mid_edge.index(x[0])+1, node_mid_edge.index(x[1])+1]
            # weight = sum((point_a-point_b)**2)
            weight = math.hypot((point_a[0] - point_b[0]),(point_a[1] - point_b[1]))

            self.mid_edges_graph.add_edge(mid_edge_idx[0], mid_edge_idx[1], weight = weight)

        return None

    def centroid(self, vertices):
        vertices.append(vertices[0])
        
        A = np.subtract([[vertices[i][0]*vertices[i+1][1]] for i in range(len(vertices)-1)] , [[vertices[i+1][0]*vertices[i][1]] for i in range(len(vertices)-1)])
        As = sum(A)/2

        x_bar = ((sum([(vertices[i+1][0] + vertices[i][0])*A[i] for i in range(len(vertices)-1)])*1/6)/As)[0]
        y_bar = ((sum([(vertices[i+1][1] + vertices[i][1])*A[i] for i in range(len(vertices)-1)])*1/6)/As)[0]
        
        vertices = vertices[: len(vertices)-1]
        return [x_bar, y_bar]


    def visual_area(self):
        for i in self.wall_edges_dict:
            wall_i = self.wall_edges_dict[i]
            points = [self.wall_graph._node[wall_i[0]]['Coordinates'], self.wall_graph._node[wall_i[1]]['Coordinates']]
            if self.wall_graph[wall_i[0]][wall_i[1]]['isWall'] == 1:
                c = 'k'
            else:
                c = '--r'
            
            plt.plot([points[0][0], points[1][0]], [points[0][1], points[1][1]], c)
                
        for i in self.area_graph._node:
            plt.text(self.area_graph._node[i]['centroid'][0], self.area_graph._node[i]['centroid'][1], str(i))
            for j in self.area_graph._node[i]['actual_vertices']:
                x = self.wall_graph._node[j]['Coordinates'][0]
                y = self.wall_graph._node[j]['Coordinates'][1]
                plt.plot(x, y, 'b*')

        plt.axis('equal')
        plt.show()

def PolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def inpolygon(xq, yq, xv, yv):
    shape = xq.shape
    xq = xq.reshape(-1)
    yq = yq.reshape(-1)
    xv = xv.reshape(-1)
    yv = yv.reshape(-1)
    q = [(xq[i], yq[i]) for i in range(xq.shape[0])]
    p = path.Path([(xv[i], yv[i]) for i in range(xv.shape[0])])
    return p.contains_points(q).reshape(shape)

if __name__ == "__main__":

    x = [ -2.868809,   -3.178729,    2.134198,    2.458878,    1.263469,    2.665492,
        -2.662194,   -2.721227,   -2.824534,   -4.093734,   -5.23011,    -9.465694,
        -10.86772,   -11.88603,   -13.8341,    -13.71604,   -11.16288,   -12.29926,
        -10.86671,    -9.907436,   -9.376143,   -9.346628,   -7.487102,   -7.516618,
        -7.35428,    -7.295245,   -8.490656,   -8.667751,   -5.273382,   -5.450479,
        -4.181279,   -4.406027,  -14.1392,    -14.18191,   -15.1707,    -18.41749,
        -18.90451,   -14.78699,   -14.38835,   -13.57665,   -14.00463,    -8.898319,
        -7.378232,   -7.153488,  -12.59924,    -7.05362,    -4.515223,   -3.290298,
        -3.231264,   -3.54458,     4.882314,    5.3103,    -17.29193,   -17.13023,
        -11.02171,   -10.98577,   -17.52548,   -17.34582,   -10.87798,   -10.78815,
        -15.8895,    -15.47627,   -10.75115,   -10.4637,    -15.6188,    -15.43913,
        -10.51638,   -10.39062,    -8.989254,   -9.204849,   -9.344091,   -9.380024,
        -9.375332,   -9.487223,   -9.541122,    7.389399,    9.186024,    8.521273,
        6.706682,   10.44366,    13.71352,    13.91114,    10.55146,     4.784107,
        10.26381,     9.976351,   13.83909,    18.35345,    13.55647,    19.26973,
        19.41238,    20.13103,    20.79578,    20.77781,    19.38692,    20.33913,
        31.92088,    32.39985,    31.03441,    31.08831,    30.01034,    20.73976,
        29.05929,    31.30507,    32.05965,    32.5088,     31.79016,    31.75423,
        28.48437,    27.94538,    26.77758,    23.93891,    24.424,      24.49379,
        22.33784,    22.07952,    17.67779,    17.28253,    17.22863,    20.44459,
        21.66629,    13.88454,    14.23571,    10.7118,     14.51698,    14.60772,
        15.09368,    11.57229,    15.98926,    12.45733,    18.13467,    18.35026,
        9.528836,    8.127469,    7.983739 ]
    y = [ 1.006395,   -3.288222,   -3.657177,    0.5193763,   0.622683,    2.851161,
        3.279147,    2.70358,     1.640994,    1.109702,    1.19825,     1.508172,
        1.626236,    1.773818,    2.009948,    4.194151,    3.928505,   -3.893306,
        -3.947384,   -3.90311,    -3.932626,   -2.707702,   -2.825767,   -4.13924,
        -0.6858368,   0.6571531,   0.7899766,  -2.707702,    0.5243301,  -0.848177,
        -0.9219675,  -4.255771,   -2.338892,   -3.422133,   -2.639953,   -3.746814,
        -10.55346,   -10.87814,    -5.638999,   -5.72755,   -10.90765,   -11.23233,
        -10.39112,    -6.767299,   -5.719471,   -6.098875,   -5.877501,   -5.715162,
        -4.696852,  -10.0242,    -10.54074,    -4.371837,    4.642451,    5.432966,
        5.289236,    6.349244,    6.870265,    8.343496,    7.858408,    9.493336,
        9.873585,   11.13122,    10.80783,    12.85598,    13.26791,    14.52555,
        14.36385,    15.2442,     15.22623,    11.63299,    10.75405,     9.963538,
        8.99235,     4.946264,    3.904222,   -2.143471,   -2.520761,   -8.126227,
        -7.874701,   -6.922489,   -7.425544,   -6.275706,   -5.808584,  -12.30466,
        -12.82568,    -8.316154,   -8.67548,    -9.056673,  -13.11704,   -13.53027,
        -11.9672,    -11.93127,   -12.03907,   -12.81161,    -9.144747,   -9.180679,
        -13.69664,    -7.392179,   -7.212516,   -6.457934,   -6.350137,   -5.613521,
        -4.477412,   -4.710973,   -4.616134,    1.67205,     1.81578,     2.426632,
        2.624261,   -4.418505,   -4.34664,    -4.113078,    3.055451,    4.555495,
        4.771091,    1.213776,    1.537168,   -2.990325,   -3.565246,   -3.87067,
        -3.942535,   -2.656886,    1.878522,    2.123971,    5.224082,    6.823626,
        12.88626,    13.19168,    23.58042,    23.99623,    23.38537,    25.64912,
        26.38573,    26.54743,    24.53521  ]

    p = np.array([x,y])
   
    s = [  1,   2,   3,   5,   1,   8,   7,  17,  16,  15,  33,  35,  36,  37,  38,  39,  39,  40,
    41,  42,  43,  44,  46,  46,  47,  48,  50,  51,  52,   2,  14,  18,  19,  13,  20,  21,
    21,  22,  23,  24,  32,  31,  30,  25,  31,  12,  11,  29,  26,  27,   5,   4,   8,   1,
    10,  12,  14,  16,  33,  18,  40,  49,  19,  21,  28,  26,  29,  32,  16,  53,  54,  55,
    56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  69,  71,  73,  68,  67,  63,
    59,  70,  72,  55,  74,   6, 124, 128, 130, 135, 133, 132, 131, 129, 127, 125, 123,  76,
    79,  78,  77,  51,  84,  85,  83,  80,  81,  82,  86,  87,  89,  90,  91,  92,  96,  88,
    92,  93,  94,  97,  98,  99, 102, 102, 100, 104, 103, 105, 106, 107, 108, 109, 113, 113,
    112, 117, 118, 119, 121, 116, 115, 115, 134, 130, 128, 126, 124,  52, 124,  77,  80,  78,
    51,  81,  88, 122, 122,  82, 120, 120, 101, 111, 110, 121, 114,  60,  17,  75]
    t = [  2,   3,   4,   1,   9,   7,   6,   7,  15,  33,  35,  36,  37,  38,  39,  34,  40,  41,
    42,  43,  44,  46,  45,  47,  48,  50,  51,  52,   3,  49,  18,  19,  13,  14,  12,  20,
    22,  28,  24,  32,  31,  30,  25,  23,  10,  11,  29,  26,  27,  28,   4,   6,   9,  10,
    11,  13,  15,  17,  34,  45,  45,  48,  20,  24,  23,  25,  30,  47,  53,  54,  55,  56,
    57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  70,  72,  74,  69,  64,  60,
    56,  71,  73,  17,  75, 124, 128, 130, 135, 134, 132, 131, 129, 127, 126, 123, 122,  79,
    78,  77,  76,  84,  85,  86,  80,  81,  82,  83,  87,  89,  90,  91,  92,  96,  95,  87,
    93,  94,  97,  98,  99, 100,  96, 101, 104, 105, 104, 106, 107, 108, 109, 110, 109, 112,
    111, 118, 119, 120, 116, 117, 116, 114, 133, 129, 127, 125, 123,  76,  77,  83,  86,  86,
    79,  87,  95, 118,  82, 102, 121, 102, 100, 110, 103, 112, 113,  73,  75,   7]
    isWall  = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]  

    con_e = np.array([s, t])

    area_edges = np.array( [[1, 2, 3, 51, 4],[6, 53, 5, 4, 51, 52, 7],
    [9, 57, 34, 56, 46, 55, 54, 5, 53, 6, 8, 58],
    [10, 59, 16, 17, 61, 60, 31, 57],[13, 14, 15, 16, 59, 11, 12],
    [18, 19, 20, 21, 22, 23, 61],
    [60, 23, 24, 68, 40, 64, 36, 63, 32],[33, 63, 35, 56],
    [37, 64, 39, 65, 38],[50, 65, 44, 66, 49],
    [66, 43, 67, 48],[47, 67, 42, 45, 55],
    [45, 41, 68, 25, 62, 30, 1, 54],[30, 62, 26, 27, 28, 29, 2],
    [70, 69, 58, 94, 71],[74, 73, 91, 75],[78, 77, 90, 79],
    [82, 81, 89, 83],[84, 89, 80, 90, 176, 93, 86, 92, 85, 88],
    [76, 91, 72, 94, 177, 95, 87, 176],
    [52, 3, 29, 158, 111, 159, 96],[28, 163, 108, 158],
    [163, 112, 113, 114, 162, 109],[110, 162, 161, 115, 160],
    [159, 160, 118, 167, 107, 157],[161, 119, 164, 116],
    [167, 168, 170, 148, 147, 166],
    [117, 164, 126, 165, 125, 133, 168],
    [120, 121, 122, 123, 124, 125, 165, 126],
    [170, 134, 171, 135, 137, 173, 172, 145, 174, 169],
    [146, 147, 148, 169, 149, 150],
    [133, 124, 127, 128, 129, 130, 131, 132, 171, 134],
    [142, 173, 137, 136, 138, 139, 140, 141],
    [144, 145, 172, 142, 143],[151, 149, 174, 144, 175, 152],
    [97, 157, 106, 156, 105, 155],[98, 155, 104, 154],
    [100, 99, 154, 103, 102, 101, 153]])

    s = [37, 36, 25, 21, 21, 24, 24, 23, 22, 26, 28, 25, 27, 30, 27, 30, 30, 30, 30, 19,  1,  2,  2,  3,
    3,  3,  3,  3,  4,  4,  4, 13,  7,  7,  9, 10, 11,  7, 18, 17, 16, 15]
    t = [38, 37, 36, 22, 25, 25, 26, 24, 23, 28, 29, 27, 28, 31, 30, 32, 34, 33, 35, 20,  2, 21,  3, 13,
    12,  8,  4, 15,  5,  7,  6, 14,  8,  9, 10, 11, 12, 13, 19, 19, 20, 20]

    con_a = np.array([s, t])

    #----------------------------------Construction Graph---------------------------------------------------------------#
    CPP = infrastructureGraph(p, con_e, isWall, True)
    CPP.initialize_graphs()
    # CPP.visual_area()