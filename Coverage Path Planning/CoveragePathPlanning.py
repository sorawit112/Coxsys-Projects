import numpy as np
from numpy import matlib
import matplotlib.pyplot as plt
from tabulate import tabulate
import networkx as nx


class infrastructureGraph:
    p = np.array([])
    con_e = np.array([])
    isWall = np.array([])
    area_edges = np.array([])
    con_a = np.array([])

    width = 0

    wall_graph = nx.Graph()
    area_graph = nx.Graph()

    wall_edges_dict = {}
    area_edges_dict = {}

    VISUAL = False

    def __init__(self, selected_points, connectivity_edges, selected_wall, selected_area_edges, connectivity_areas, wall_offset_width, visual = False):
        self.p = selected_points
        self.con_e = connectivity_edges
        self.isWall = selected_wall
        self.area_edges = selected_area_edges
        self.con_a = connectivity_areas
        self.width = wall_offset_width
        self.VISUAL = visual

    def initialize_graphs(self):
        self.wallGraphInit()
        self.areaGraphInit()

        return None

    def wallGraphInit(self):
        name_Nodes = range(1, self.p.shape[1]+1)
        name_Edges = range(1, self.con_e[0].shape[0]+1)

        # EdgeTable = tabulate([[np.transpose(self.con_e), np.transpose(name_Edges), np.transpose(self.isWall)]], headers=["EndNodes", "Number", "isWall"], tablefmt="fancy_grid")
        # NodeTable = tabulate([[np.transpose(name_Nodes), np.transpose(p)]], headers=["Number", "Coordinates"], tablefmt="fancy_grid")
        # print(NodeTable)
        # print(EdgeTable)
        
        # create nodes
        for node in name_Nodes:
            self.wall_graph.add_node(node, Coordinates = [self.p[0][node-1], self.p[1][node-1]])

        # create edges
        for edge in name_Edges:
            self.wall_edges_dict[edge] = [self.con_e[0][edge-1], self.con_e[1][edge-1]]
            self.wall_graph.add_edge(self.con_e[0][edge-1], self.con_e[1][edge-1], Number = edge, isWall = self.isWall[0][edge-1])

        return None
    
    def areaGraphInit(self):
        num_area = self.area_edges.shape[0]
        vertices = self.detect_vertices(num_area)

        name_Nodes = range(1, num_area+1)
        num_edges = self.con_a.shape[1]
        # NodeTable = tabulate([np.transpose(name_Nodes), np.transpose(centroid), np.transpose(area), np.transpose(vertices), np.transpose(inner_vertices), np.transpose(corner_vertices), np.transpose(area_edges), np.transpose(theta), np.transpose(alpha), np.transpose(ver_opp), np.transpose(rho_offset)],headers=["Number","Centroid","Area","Vertices","InnerVertices","Corner","Edges","Theta","Alpha","VerOpp","RhoOffset"], tablefmt="fancy_grid")
        
        # create nodes
        for node in name_Nodes:
            self.area_graph.add_node(node, centroid = [0, 0], area = 0, vertices = vertices[node-1], inner_vertices = [], corner_vertices = [], area_edges = self.area_edges[node-1], theta = [], alpha = [], ver_opp = [], rho_offset = [])

        # create edges
        for edge in range(num_edges):
            self.area_edges_dict[edge] = [self.con_a[0][edge], self.con_a[1][edge]]
            self.area_graph.add_edge(self.con_a[0][edge], self.con_a[1][edge], distance = 0)


        return None

    def areaNodes_properties(self):
        
        # calculate offset
        for i in self.area_graph._node:
            coordinates = [self.wall_graph._node[ver]['Coordinates'] for ver in self.area_graph._node[i]['vertices']]
            self.area_graph._node[i]['area'] = PolyArea([cor[0] for cor in coordinates], [cor[1] for cor in coordinates])
            num_ver = len(self.area_graph._node[i]['vertices'])
            centroid_i = list(np.sum(coordinates, axis=0)/num_ver)
            self.area_graph._node[i]['centroid'] = centroid_i

            ver_shift = np.subtract(coordinates, centroid_i)
            theta_i = np.zeros((num_ver))
            alpha_i = np.zeros((num_ver))
            rho = np.zeros((num_ver))
            rho_offset_i = np.zeros((num_ver))
            rho_temp = np.zeros((num_ver))
            p = np.zeros((num_ver, 2))
            p_temp = np.zeros((num_ver, 2))

            for j in range(num_ver):
                point_a = ver_shift[j]
                if j+1 != num_ver:
                    point_b = ver_shift[j+1]
                else:
                    point_b = ver_shift[0]

                theta_i[j] = np.mod(np.arctan2(point_b[1] - point_a[1], point_b[0] - point_a[0]), 2*np.pi)
                rho[j] = -np.cos(theta_i[j])*point_b[1] + np.sin(theta_i[j])*point_b[0]
                rho_offset_i[j] = rho[j] - self.width*0.2
                rho_temp[j] = rho[j] - self.width/(np.sqrt(2) + 0.1)
            
            # print(theta_i)
            # print(rho)
            # print(rho_offset_i)
            # print(rho_temp)

            for j in range(num_ver):
                theta_a = theta_i[j]
                rho_a = rho_offset_i[j]
                rho_c = rho_temp[j]

                if j+1 != num_ver:
                    idx = j+1
                else:
                    idx = 0

                theta_b = theta_i[idx]
                rho_b = rho_offset_i[idx]
                rho_d = rho_temp[idx]
                theta_b = theta_b + (theta_b < theta_a)*2*np.pi
                alpha_i[idx] = (theta_a + theta_b)/2 + np.pi/2
                p_equ = np.matmul(np.linalg.inv([[np.sin(theta_a), -np.cos(theta_a)],[np.sin(theta_b), -np.cos(theta_b)]]), [[rho_a], [rho_b]])

                p[idx][0] = p_equ[0]
                p[idx][1] = p_equ[1]

                p_temp_equ = np.matmul(np.linalg.inv([[np.sin(theta_a), -np.cos(theta_a)],[np.sin(theta_b), -np.cos(theta_b)]]), [[rho_c], [rho_d]])

                p_temp[idx][0] = p_temp_equ[0]
                p_temp[idx][1] = p_temp_equ[1]

            self.area_graph._node[i]['theta'] = theta_i
            self.area_graph._node[i]['alpha'] = alpha_i
            self.area_graph._node[i]['rho_offset'] = rho_offset_i
            self.area_graph._node[i]['inner_vertices'] = np.transpose(np.add(p, self.area_graph._node[i]['centroid']))
            self.area_graph._node[i]['corner_vertices'] = np.transpose(np.add(p_temp, self.area_graph._node[i]['centroid']))

            # print(theta_i)
            # print(alpha_i)
            # print(rho_offset_i)
            # print(self.area_graph._node[i]['inner_vertices'])
            # print(self.area_graph._node[i]['corner'])

        # calculate distance btw connectivity area 
        num_edges = len(self.area_edges)
        for i in range(num_edges):
            cm_ab = [self.area_graph._node[self.area_edges_dict[i][0]]['centroid'], self.area_graph._node[self.area_edges_dict[i][1]]['centroid']]
            self.area_graph[self.area_edges_dict[i][0]][self.area_edges_dict[i][1]]['distance'] = np.sqrt((cm_ab[0][0]-cm_ab[1][0])**2 + (cm_ab[0][1]-cm_ab[1][1])**2)
        
        # find opposite vertices
        for i in self.area_graph._node:
            ver_list = self.area_graph._node[i]['vertices']
            alpha_i = self.area_graph._node[i]['alpha']

            num_ver = len(ver_list)
            ver_opp_i = np.zeros(num_ver)

            for j in range(num_ver):
                ver_start = ver_list[j]
                wall_sub = [self.wall_graph._node[x]['Coordinates'] for x in ver_list]

                c = np.cos(-alpha_i[j])
                s = np.sin(-alpha_i[j])
                R = [[c, -s], [s, c]]

                offset = [self.wall_graph._node[ver_list[j]]['Coordinates']]
                new_area = np.matmul(R, (np.transpose(wall_sub)-np.transpose(offset)))
                
                max_idx = [idx for idx,val in enumerate(new_area[0]) if val == max(new_area[0])][0]
                ver_opp_i[j] = ver_list[max_idx]
            
            self.area_graph._node[i]['ver_opp'] = ver_opp_i

        return None

    def detect_vertices(self, num_area):
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

        return vertices

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

            for i in self.area_graph._node:
                plt.text(self.area_graph._node[i]['centroid'][0], self.area_graph._node[i]['centroid'][1], str(i))
                num_ver = len(self.area_graph._node[i]['vertices'])
                corner_i = self.area_graph._node[i]['corner_vertices']
                alpha_i = self.area_graph._node[i]['alpha']

            #     for j in range(num_ver):
            #         plt.plot(corner_i[0][j], corner_i[1][j], 'mo')
            #         L = 1
            #         plt.plot([corner_i[0][j], corner_i[0][j]+L*np.cos(alpha_i[j])], [corner_i[1][j], corner_i[1][j]+L*np.sin(alpha_i[j])], 'b')

            plt.axis('equal')
            plt.show()
        
class coveragePathPlanning(infrastructureGraph):
    start_area = 7
    goal_area = 1
    cleanning_area = [2, 3, 5, 6, 7, 8]

    def __init__(self, selected_points, connectivity_edges, selected_wall, selected_area_edges, connectivity_areas, wall_offset_width, visual = False):
        super().__init__(selected_points, connectivity_edges, selected_wall, selected_area_edges, connectivity_areas, wall_offset_width, visual)

    def generate_robotWorkFlowPath(self, start_area, goal_area, cleaning_area):
        
        pass

    def generate_viapoints(self):
        pass
    
    def generate_areaPath(self):
        k = 1
        current_graph = [self.start_area]
        candidate = []
        max_iter = 5
        break_loop = False

        while (~break_loop):
            while (k <= max_iter):
                temp_graph = []
                temp_edges = []
                
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
                break_loop = True

        # find minimum path
        dist_path = []
        for candi_i in candidate:
            dist_sum = 0
            for i in range(len(candi_i)-2):
                dist_sum = dist_sum + self.area_graph[candi_i[i]][candi_i[i+1]]




def PolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

if __name__ == "__main__":
#----------------data from user interface--------------------------------------------------------------------------#
 
# vertices (selected points)
    x = np.array([-8, -4, 0, 2,  0,  -5,  5, 10,  3,   4,  6, 11, 17,  13, 21])
    y = np.array([ 5, 0, 9, 5, -4, -10, -1,  2, -8, -15, -6, -8,  4, -15, -5])
    p = np.array([x, y])

# connectivity edges (drag and draw)
    s = np.array([1, 1, 2, 2, 2, 3, 4, 4, 5, 5, 5,  6, 7,  7,  8,  8,  9,  9, 10, 10, 11, 12, 12, 13])
    t = np.array([2, 3, 4 ,5 ,6 ,4, 7, 8, 6, 7, 9, 10, 8, 11, 12, 13, 10, 11, 12, 14, 12, 14, 15, 15])
    con_e = np.array([s, t])
    isWall = np.array([[1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1]])

# selected_areas (consecutive_edges only)
    area_edges = np.array([[1, 3, 6, 2], [3, 4, 10, 7], [4, 5, 9], [7, 13, 8], [9, 12, 17, 11], [13, 14, 21, 15], [17, 19, 21, 18], [19, 20, 22], [15, 23, 24, 16]])

# connectivity areas (drag and draw)
    s = np.array([1, 2, 2, 3, 4, 5, 6, 6, 7])
    t = np.array([2, 3, 4, 5, 6, 7, 7, 9, 8])
    con_a = np.array([s, t])

#----------------------------------Construction Graph---------------------------------------------------------------#
    CPP = coveragePathPlanning(p, con_e, isWall, area_edges, con_a, 0.75, True)
    CPP.initialize_graphs()
    CPP.areaNodes_properties()
    # CPP.visualize_area()
    CPP.generate_areaPath()
    
    print("")