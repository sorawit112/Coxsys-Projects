import numpy as np
from numpy import matlib
import matplotlib.pyplot as plt
from tabulate import tabulate
import networkx as nx


class infrastructureGraph:
    energy = 67   # to kill [J/m^2]

    width_inner = 0.3 
    power_inner = 0.76

    width_outer = 1.3  
    power_outer = 0.53

    offset_ratio = 1.5

    grid_width = (width_outer-width_inner)
    edge_width = (width_outer-width_inner)
    power_actual = min(power_inner,power_outer)

    duration = energy/power_actual
    minute = duration/60


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
    mid_edges_dict = {}

    VISUAL = False

    def __init__(self, selected_points, connectivity_edges, selected_wall, selected_area_edges, connectivity_areas, visual = False):
        self.p = selected_points
        self.con_e = connectivity_edges
        self.isWall = selected_wall
        self.area_edges = selected_area_edges
        self.con_a = connectivity_areas
        self.VISUAL = visual

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
            self.wall_edges_dict[edge] = [self.con_e[0][edge-1], self.con_e[1][edge-1]]
            self.wall_graph.add_edge(self.con_e[0][edge-1], self.con_e[1][edge-1], Number = edge, isWall = self.isWall[edge-1], midEdge = np.add([self.wall_graph._node[self.con_e[0][edge-1]]['Coordinates'][0], self.wall_graph._node[self.con_e[0][edge-1]]['Coordinates'][1]],[self.wall_graph._node[self.con_e[1][edge-1]]['Coordinates'][0], self.wall_graph._node[self.con_e[1][edge-1]]['Coordinates'][1]])/2)

        return None
    
    def areaGraphInit(self):
        num_area = self.area_edges.shape[0]
        vertices = []

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
        # NodeTable = tabulate([np.transpose(name_Nodes), np.transpose(centroid), np.transpose(area), np.transpose(vertices), np.transpose(inner_vertices), np.transpose(corner_vertices), np.transpose(area_edges), np.transpose(theta), np.transpose(alpha), np.transpose(ver_opp), np.transpose(rho_offset)],headers=["Number","Centroid","Area","Vertices","InnerVertices","Corner","Edges","Theta","Alpha","VerOpp","RhoOffset"], tablefmt="fancy_grid")
        
        # create nodes
        for node in name_Nodes:
            self.area_graph.add_node(node, centroid = [0, 0], vertices = vertices[node-1], actual_vertices = [], area_edges = self.area_edges[node-1])

        # create edges
        for edge in range(num_edges):
            self.area_edges_dict[edge] = [self.con_a[0][edge], self.con_a[1][edge]]
            self.area_graph.add_edge(self.con_a[0][edge], self.con_a[1][edge])
    
        # calculate nodes properties
        caltimes = 0
        while caltimes < 2: 
            for i in self.area_graph._node:
                coordinates = [self.wall_graph._node[ver]['Coordinates'] for ver in self.area_graph._node[i]['vertices']]
                num_ver = len(self.area_graph._node[i]['vertices'])
                centroid_i = list(np.sum(coordinates, axis=0)/num_ver)
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

                num_actual_ver = len(actual_theta)

                self.area_graph._node[i]['actual_vertices'] = actual_vertices

                actual_coordinates = [self.wall_graph._node[x]['Coordinates'] for x in actual_vertices]
                num_actual_ver = len(actual_vertices)
                self.area_graph._node[i]['centroid'] = sum(np.array(actual_coordinates))/num_actual_ver

            caltimes = caltimes+1
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
            weight = sum((point_a-point_b)**2)

            self.mid_edges_graph.add_edge(mid_edge_idx[0], mid_edge_idx[1], weight = weight)

        return None

    def centroid(self, vertices_i):
        vertices = [vertices_i , vertices_i]
        A = vertices


    def visualize_area(self):
        if self.VISUAL:
            for i in self.wall_edges_dict:
                wall_i = self.wall_edges_dict[i]
                points = [self.wall_graph._node[wall_i[0]]['Coordinates'], self.wall_graph._node[wall_i[1]]['Coordinates']]
                if self.wall_graph[wall_i[0]][wall_i[1]]['isWall'] == 1:
                    c = 'k'
                else:
                    c = 'r'
                
                plt.plot([points[0][0], points[1][0]], [points[0][1], points[1][1]], c)

                point_mid = self.wall_graph[wall_i[0]][wall_i[1]]['midEdge']
                plt.plot(point_mid[0], point_mid[1], 'k.')
                plt.text(point_mid[0], point_mid[1], str(i))
                

            for i in self.area_graph._node:
                plt.text(self.area_graph._node[i]['centroid'][0], self.area_graph._node[i]['centroid'][1], str(i))
                for j in self.area_graph._node[i]['actual_vertices']:
                    x = self.wall_graph._node[j]['Coordinates'][0]
                    y = self.wall_graph._node[j]['Coordinates'][1]
                    plt.plot(x, y, 'b*')
                
            plt.axis('equal')
            plt.show()
        
class coveragePathPlanning(infrastructureGraph):
    start_area = 1
    goal_area = 7
    cleanning_area = [2, 3, 4, 5, 8, 9]
    area_path = []
    clean_path = []
    mid_edge_path_idx = []

    def __init__(self, selected_points, connectivity_edges, selected_wall, selected_area_edges, connectivity_areas, visual = False):
        super().__init__(selected_points, connectivity_edges, selected_wall, selected_area_edges, connectivity_areas, visual)

    def generate_robotWorkFlowPath(self, start_area, goal_area, cleaning_area):
        
        pass

    def generate_viapoints(self):
        pass
    
    def generate_areaPath(self):
        k = 1
        current_graph = [self.start_area]
        candidate = []
        max_iter = 4

        while True:
            while (k <= max_iter):
                temp_graph = []
                
                if type(current_graph[0]) == int:
                    neighbor_list = [[n] for n in self.area_graph.neighbors(current_graph[0])]
                    temp_graph = np.concatenate((matlib.repmat(current_graph, len(neighbor_list),1), neighbor_list), axis = 1)
                else: 
                    graph_i = []
                    num = len(current_graph[0])
                    for i in range(len(current_graph)):
                        neighbor_list = [[n] for n in self.area_graph.neighbors(current_graph[i][num-1])]
                        graph_i = np.concatenate((matlib.repmat(current_graph[i], len(neighbor_list), 1), neighbor_list),axis = 1)

                        if i == 0:
                            temp_graph = graph_i
                        else:
                            temp_graph = np.concatenate((temp_graph, graph_i))

                current_graph = temp_graph
                num = len(current_graph[0])
                candidate = [path_list for path_list in current_graph if (path_list[num-1] == self.goal_area) and (set(self.cleanning_area).issubset(set(path_list)))]
                
                k = k+1

            if len(candidate) == 0:
                max_iter = max_iter+1
            else:
                break

        # find minimum path
        dist_path = []
        print(candidate)
        for candi_i in candidate:
            path_edge = [list(set([e for e in self.area_graph._node[candi_i[x]]['area_edges']]).intersection(set([e for e in self.area_graph._node[candi_i[x+1]]['area_edges']])))[0] for x in range(len(candi_i)-1)]
            node_idx = [[self.mid_edges_graph._node[y]['edge_Number'] for y in self.mid_edges_graph._node].index(x)+1 for x in path_edge]
            dist_temp = 0

            for idx in range(len(node_idx)-1):
                current_node = node_idx[idx]
                next_node = node_idx[idx+1]
                if next_node - current_node != 0:
                    dist_temp = dist_temp + self.mid_edges_graph[current_node][next_node]['weight']
            dist_path.append(dist_temp)

        min_dist = min(dist_path)
        min_dist_idx = dist_path.index(min_dist)

        path = candidate[min_dist_idx]

        cleanList = self.cleanning_area
        path_clean = []
        for i in range(len(path)):
            if path[i] in cleanList:
                path_clean.append(1)
                cleanList.pop(cleanList.index(path[i]))
            else:
                path_clean.append(0)

        mid_edge_path = [list(set([e for e in self.area_graph._node[path[x]]['area_edges']]).intersection(set([e for e in self.area_graph._node[path[x+1]]['area_edges']])))[0] for x in range(len(path)-1)]
        
        self.area_path = path
        self.clean_path = path_clean
        self.mid_edge_path_idx = [[self.mid_edges_graph._node[y]['edge_Number'] for y in self.mid_edges_graph._node].index(x)+1 for x in mid_edge_path]
               
        return None

if __name__ == "__main__":
#----------------data from user interface--------------------------------------------------------------------------#
 
# vertices (selected points)
    x = [-75, -75, -75+110, -75+110, -75+110, -75+110, -75+110,  -75+110+45, -75+110+45, -75+110+45, -75+110+45, -75+110+50, -75+110+50, -75+110+50, -75+100+45+22, -75+100+45+22, -75+110+50+25, -75+110+50+25+45, -75+110+50+80, -75+110+192-40, -75+110+192-40, -75+110+192-40, -75+110+192, -75+110+192, -75+110+192, -75+110+192, -75+110+192+550, -75+110+192+550]
    y = [10-40, 10, 10-40, 10, 10-40+100, 10-40+100+50, 10-40+100+50+562, 10-40+100+50, 10-40+100+50+45, 10-40+100+50+45+50, 10-40+100+50+562, 10-40, 10-40+60, 10-40+100, 10-40+100+50+45, 10-40+100+50+45+50, 10-40, 10-40-45, 10-40+60-50, 10-40-45, 10-40+60-50, 10-40+100, 10-40-45, 10-40, 10-40+100, 10-40+150, 10-40-45, 10-40]
    # x = np.array([-8, -4, 0, 2,  0,  -5,  5, 10,  3,   4,  6, 11, 17,  13, 21])
    # y = np.array([ 5, 0, 9, 5, -4, -10, -1,  2, -8, -15, -6, -8,  4, -15, -5])
    p = np.array([x, y])*4/90

# connectivity edges (drag and draw)
    s =      [1, 1, 2, 3, 3, 4, 5, 5,  6, 6, 7,  8, 8,  9,  9,  10, 10, 12, 12, 13, 13, 14, 15, 17, 17, 18, 19, 20, 20, 21, 22, 23, 23, 24, 24, 25, 27]
    t =      [2, 3, 4, 4,12, 5, 6, 14, 7, 8, 11, 9, 26, 10, 15, 11, 16, 13, 17, 14, 19, 22, 16, 18, 19, 20, 21, 21, 23, 22, 25, 24, 27, 25, 28, 26, 28]
    isWall = [1, 1, 1, 0, 1, 1, 1, 0,  1, 0, 1,  1, 1,  0,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,  0,  1,  1,  0,  1,  1,  0,  0,  1,  1,  1,  1,  1]

    # s = np.array([1, 1, 2, 2, 2, 3, 4, 4, 5, 5, 5,  6, 7,  7,  8,  8,  9,  9, 10, 10, 11, 12, 12, 13])
    # t = np.array([2, 3, 4 ,5 ,6 ,4, 7, 8, 6, 7, 9, 10, 8, 11, 12, 13, 10, 11, 12, 14, 12, 14, 15, 15])
    # isWall = np.array([1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1])

    con_e = np.array([s, t])

# selected_areas (consecutive_edges only)
    area_edges = np.array([[2, 4, 3, 1],[5, 18, 20, 8, 6, 4],[8, 22, 31, 36, 13, 10, 7],[10, 12, 14, 16, 11, 9],[15, 23, 17, 14],[19, 25, 21, 18],[24, 26, 28, 27, 25],[29, 32, 34, 31, 30, 28],[33, 37, 35, 32]])
    # area_edges = np.array([[1, 3, 6, 2], [3, 4, 10, 7], [4, 5, 9], [7, 13, 8], [9, 12, 17, 11], [13, 14, 21, 15], [17, 19, 21, 18], [19, 20, 22], [15, 23, 24, 16]])
    
# connectivity areas (drag and draw)
    s = [1, 2, 2, 3, 3, 4, 6, 7, 8]
    t = [2, 3, 6, 4, 8, 5, 7, 8, 9]

    # s = [1, 2, 2, 3, 4, 5, 6, 6, 7]
    # t = [2, 3, 4, 5, 6, 7, 7, 9, 8]
    con_a = np.array([s, t])

#----------------------------------Construction Graph---------------------------------------------------------------#
    CPP = coveragePathPlanning(p, con_e, isWall, area_edges, con_a, True)
    CPP.initialize_graphs()
    CPP.generate_areaPath()
    CPP.visualize_area()
    
    print("")