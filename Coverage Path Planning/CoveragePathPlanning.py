from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp

from matplotlib import animation
from matplotlib import path
from numpy import matlib

from scipy.cluster.vq import kmeans,vq
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

import datetime
import math

##############################################################################################################################
#####                                                                                                                    #####
##### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         Construct Graph Unit        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #####
#####                                                                                                                    #####
##############################################################################################################################

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

    def __init__(self, selected_points, connectivity_edges, selected_wall, selected_area_edges, connectivity_areas, visual = False):
        self.p = selected_points
        self.con_e = connectivity_edges
        self.isWall = selected_wall
        self.area_edges = selected_area_edges
        self.con_a = connectivity_areas
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
        vertices = []

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


##############################################################################################################################
#####                                                                                                                    #####
##### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         Coverage path planning Unit         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #####
#####                                                                                                                    #####
##############################################################################################################################

class coveragePathPlanning(infrastructureGraph):
    max_cleanning_radius = 1.5
    max_iteration_area_search = 5
    energy_to_kill = 67
    coef = [[0.1431328927561467724],
          [-11.0300107335172262],
              [80.2553874972875],
       [-209.437642843722228422],
            [197.34549417652961],
          [-3.61325901328181954]]

    circle_shift = 30 #degree
    animate_times = 400 #ms    

    def __init__(self, selected_points, connectivity_edges, selected_wall, selected_area_edges, connectivity_areas, visual = False):        
        super().__init__(selected_points, connectivity_edges, selected_wall, selected_area_edges, connectivity_areas, visual)

        self.time_kill = []
        self.start_area = []
        self.goal_area = []

        self.area_path = []
        self.clean_path = []
        self.cleanning_area = []
        self.mid_edge_path_idx = []

        self.robot_path = []
        self.robot_yaw = []
        self.robot_duration = []
        self.robot_yaw = []
        self.robot_qw = []
        self.robot_clean_radius = []
        self.robot_clean_status = []

    ##############################################################################################################################
    #####                                                                                                                    #####
    #####                                             main function                                                          #####
    #####                                                                                                                    #####
    ##############################################################################################################################

    def generate_robotWorkFlowPath(self, start_pose, goal_pose, cleaning_area):        

        self.fleet_management( start_pose, goal_pose, cleaning_area)

        ## generate via-points
        n_robot = len(start_pose)
        flag_shift = False
        for r in range(n_robot):
            num_node_path = len(self.area_path[r])
            current_point = list(start_pose[r])
            self.time_kill.append([])
            self.robot_path.append([current_point])
            self.robot_duration.append([0])
            self.robot_yaw.append([0])
            self.robot_qw.append([yaw_to_qw(0)])
            self.robot_clean_status.append([0])
            self.robot_clean_radius.append([0])

            time = 0

            for i in range(num_node_path):
                
                if self.clean_path[r][i]:
                    if i != num_node_path-1:
                        if self.clean_path[r][i+1] == 1 or i+1 == len(self.mid_edge_path_idx[r]):
                            edge = self.wall_edges_dict[self.mid_edges_graph._node[self.mid_edge_path_idx[r][i]]['edge_Number']]
                            end_point = self.wall_graph[edge[0]][edge[1]]['midEdge']
                        else:
                            edge = self.wall_edges_dict[self.mid_edges_graph._node[self.mid_edge_path_idx[r][i+1]]['edge_Number']]
                            end_point = self.wall_graph[edge[0]][edge[1]]['midEdge']

                    else:
                        end_point = self.area_graph._node[self.area_path[r][i]]['centroid']

                    [via_points, region_radius] = self.viaPointGenerator([self.wall_graph._node[ver]['Coordinates'] for ver in self.area_graph._node[self.area_path[r][i]]['actual_vertices']], current_point, end_point)

                    time = time + self.calculate_time_kill(region_radius)

                    for i in range(len(via_points)):
                        delta_yaw = np.arctan2(via_points[i][1]-self.robot_path[r][len(self.robot_path[r])-1][1], via_points[i][0]-self.robot_path[r][len(self.robot_path[r])-1][0])
                        self.robot_yaw[r].append(delta_yaw)
                        self.robot_qw[r].append(yaw_to_qw(delta_yaw))

                        self.robot_path[r].append(via_points[i])
                        self.robot_duration[r].append(self.calculate_time_kill(region_radius[i]))

                        self.robot_clean_radius[r].append(region_radius[i])
                        self.robot_clean_status[r].append(1)
                        
                    current_point = via_points[len(via_points)-1]

                else:
                    if i != num_node_path-1:
                        if i != 0 and self.clean_path[r][i-1]:
                            edge1 = self.wall_edges_dict[self.mid_edges_graph._node[self.mid_edge_path_idx[r][i-1]]['edge_Number']]
                            edge2 = self.wall_edges_dict[self.mid_edges_graph._node[self.mid_edge_path_idx[r][i]]['edge_Number']]
                            next_point = [self.wall_graph[edge1[0]][edge1[1]]['midEdge'], self.wall_graph[edge2[0]][edge2[1]]['midEdge']]
                        
                        else:
                            edge = self.wall_edges_dict[self.mid_edges_graph._node[self.mid_edge_path_idx[r][i]]['edge_Number']]
                            next_point = [self.wall_graph[edge[0]][edge[1]]['midEdge']]
                    else:
                        next_point = [self.area_graph._node[self.area_path[r][i]]['centroid']]
                        
                    for n in next_point:
                        delta_yaw = np.arctan2(n[1]-self.robot_path[r][len(self.robot_path[r])-1][1], n[0]-self.robot_path[r][len(self.robot_path[r])-1][0])
                        self.robot_yaw[r].append(delta_yaw)
                        self.robot_qw[r].append(yaw_to_qw(delta_yaw))
                        
                        self.robot_path[r].append(list(n))
                        self.robot_duration[r].append(0)
                        self.robot_clean_radius[r].append(0)
                        self.robot_clean_status[r].append(0)

                    current_point = next_point[-1]
            
            self.time_kill[r] = time

            hours = np.floor(time/3600)
            minutes = np.floor((time%3600)/60)
            seconds = np.floor(time - 3600*hours - 60*minutes)
            area = sum([self.area_graph._node[area]['area'] for area in self.area_path[r]])
            
            performance = time/area

            print('---------------------- robot number ' + str(r+1)+ ' ------------------------------\n')
            print('robot path : ' + str(len(self.robot_path[r])) + ' points\n')
            for i in range(len(self.robot_path[r])):
                print('via-point '+ str(i) +' : ' + str(self.robot_path[r][i]) + ' yaw : ' + str(self.robot_yaw[r][i]) +'\n' )

            print('Time Spent : ' + str(hours) + ' Hours ' + str(minutes) + ' Minutes ' + str(seconds) + ' Seconds\n')
            print('Performance : ' + str(round(performance,2)) + ' s/m^2\n')
            print('---------------------- done robot number ' + str(r+1)+ ' -------------------------\n')
        
        self.writeText()

        if self.VISUAL:
            self.visualize_path()
    
    ##############################################################################################################################
    #####                                                                                                                    #####
    #####                                          coverage path planning part                                               #####
    #####                                                                                                                    #####
    ##############################################################################################################################

    def fleet_management(self, start_pose, goal_pose, cleaning_area):
        n_robot = len(start_pose)
        
        ## clustering group of zone with num of robots
        # generate data
        self.K_data = np.vstack([self.area_graph._node[a]['centroid'] for a in self.area_graph._node])
        
        # compute K-Means with K = n_robot
        self.K_centroids,_ = kmeans(self.K_data, n_robot)
        self.K_idx,_ = vq(self.K_data, self.K_centroids) 

        # assight K-Means to each robot
        self.K_robot = []
        selected_centroid = []
        for r in range(n_robot):
            dist = []
            for c in self.K_centroids:
                dist.append(int(math.hypot((start_pose[r][0] - c[0]),(start_pose[r][1] - c[1]))))
            
            sorted_dist = sorted(dist)
            for i in sorted_dist:
                if dist.index(i) not in selected_centroid:
                    self.K_robot.append(dist.index(i))
                    selected_centroid.append(dist.index(i))
                    break

        # manage cleaning area
        for r in range(n_robot):
            start_point = start_pose[r]
            goal_point = goal_pose[r]

            # detect area for start pose and goal pose
            i = 1
            flag_start = False
            flag_goal = False

            while i <= len(self.area_graph._node):
                ver_list = self.area_graph._node[i]['actual_vertices']
                coor = np.transpose([self.wall_graph._node[x]['Coordinates'] for x in ver_list])
                if inpolygon(np.array([start_point[0]]), np.array([start_point[1]]), coor[0], coor[1])[0] and not flag_start:
                    self.start_area.append(i)
                    flag_start = True
                
                if inpolygon(np.array([goal_point[0]]), np.array([goal_point[1]]), coor[0], coor[1])[0] and not flag_goal:
                    self.goal_area.append(i)
                    flag_goal = True

                if flag_start and flag_goal:
                    break

                i = i + 1
            
            self.cleanning_area.append([])

        # assign responsible area from devided area
        for a in cleaning_area:
            cluster = self.K_idx[a-1]

            index_robot = self.K_robot.index(cluster)
            
            self.cleanning_area[index_robot].append(a)
        
        # share load
        if n_robot > 1:
            flag_balance = False

            average_load = sum([self.area_graph._node[a]['area'] for a in cleaning_area])/len(cleaning_area)

            while not flag_balance:
                unequal_load_list = []
                for r in range(n_robot):   # calculate unequal load of every robots
                    if len(self.cleanning_area[r]) != 0:
                        unequal_load_list.append(sum([self.area_graph._node[a]['area'] for a in self.cleanning_area[r]]))
                    else:
                        unequal_load_list.append(0)
                max_unequal_index = unequal_load_list.index(max(unequal_load_list)) 

                share_load_list = []
                current_unequal_load = self.sum_unequal_load(cleaning_area)
                
                for i in range(len(self.cleanning_area[max_unequal_index])):
                    area = self.cleanning_area[max_unequal_index][i]
                    self.cleanning_area[max_unequal_index].pop(i)     
                    share_load_list.append([])    
                    for r in range(n_robot):
                        if r != max_unequal_index:
                            self.cleanning_area[r].append(area)
                            share_load_list[len(share_load_list)-1].append(self.sum_unequal_load(cleaning_area))
                            self.cleanning_area[r].pop(len(self.cleanning_area[r])-1)
                        else:
                            share_load_list[len(share_load_list)-1].append(10000000)
                    self.cleanning_area[max_unequal_index].insert(i,area)

                min_share_load_robot = min(share_load_list)
                robot_index = min(share_load_list).index(min(min_share_load_robot))

                min_share_load_area = min_share_load_robot[robot_index]
                area_index = list(np.transpose(share_load_list)[robot_index]).index(min_share_load_area)

                # find closet area to share
                distance_weight_list = []
                for i in self.cleanning_area[max_unequal_index]:
                    from_node = self.area_graph._node[i]['centroid']
                    dist = 0
                    if len(self.cleanning_area[robot_index]) != 0:
                        for j in self.cleanning_area[robot_index]:
                            path, dist_ = self.find_shortest_path(i, j)
                            dist += dist_
                    else:
                        path, dist_ = self.find_shortest_path(i, self.start_area[robot_index])
                        dist += dist_
                    distance_weight_list.append(dist)
                
                closet_index = distance_weight_list.index(min(distance_weight_list))
                self.cleanning_area[robot_index].append(self.cleanning_area[max_unequal_index][closet_index])
                self.cleanning_area[max_unequal_index].pop(closet_index)
                
                next_unequal_load = self.sum_unequal_load(cleaning_area)
                if next_unequal_load - current_unequal_load >= 0:   #if share sum_load is increse then break
                    self.cleanning_area[max_unequal_index].append(self.cleanning_area[robot_index][-1])
                    self.cleanning_area[robot_index].pop(len(self.cleanning_area[robot_index])-1)
                    break

        self.generate_areaPath()
    
    def sum_unequal_load(self, cleaning_area):
        n_robot = len(self.start_area)

        unequal_load = 0        
        average_load = sum([self.area_graph._node[a]['area'] for a in cleaning_area])/len(cleaning_area)
        
        for r1 in range(n_robot):
            for r2 in range(r1, n_robot):
                if len(self.cleanning_area[r1]) != 0:
                    load_r1 = sum([self.area_graph._node[a]['area'] for a in self.cleanning_area[r1]])
                else:
                    load_r1 = 0
                
                if len(self.cleanning_area[r2]) != 0:
                    load_r2 = sum([self.area_graph._node[a]['area'] for a in self.cleanning_area[r2]])
                else:
                    load_r2 = 0
                
                unequal_load += abs(load_r1 - load_r2)

        return unequal_load

    def generate_areaPath(self):
        n_robot = len(self.start_area)

        for r in range(n_robot):
            self.area_path.append([])
            self.clean_path.append([])
            self.mid_edge_path_idx.append([])

            cleaning_area = []
            for i in range(len(self.cleanning_area[r])):
                cleaning_area.append(self.cleanning_area[r][i])

            data = {}
            data['location'] = []
            locations = []
            if self.start_area[r] not in cleaning_area:
                data['location'] += [tuple(self.area_graph._node[self.start_area[r]]['centroid'])]
                locations.append(self.start_area[r])
            else:
                data['location'] += [tuple(self.area_graph._node[self.start_area[r]]['centroid'])]
                locations.append(self.start_area[r])
                cleaning_area.pop(cleaning_area.index(self.start_area[r]))

            for i in cleaning_area:
                data['location'] += [tuple(self.area_graph._node[i]['centroid'])]
                locations.append(i)
            
            if self.goal_area[r] not in cleaning_area:
                data['location'] += [tuple(self.area_graph._node[self.goal_area[r]]['centroid'])]
                locations.append(self.goal_area[r])
            else:
                data['location'] += [tuple(self.area_graph._node[self.goal_area[r]]['centroid'])]
                locations.append(self.goal_area[r])
                cleaning_area.pop(cleaning_area.index(self.goal_area[r]))

            data['location'] += [tuple([10000,10000])]
            locations.append('d')

            data['num_vehicles'] = 1
            data['depot'] = 1

            distances = {}
            for from_counter, from_node in enumerate(locations):
                distances[from_counter] = {}
                for to_counter, to_node in enumerate(locations):
                    if from_counter == to_counter:
                        distances[from_counter][to_counter] = 0
                    else:
                        if from_counter != len(locations)-1:
                            if to_counter != len(locations)-1:
                                if from_node == to_node:
                                    distances[from_counter][to_counter] = 100000000
                                else:    
                                    path,dist_ = self.find_shortest_path(from_node, to_node)
                                    distances[from_counter][to_counter] = dist_ 
                            else:
                                if from_counter == 0 or from_counter == len(locations)-2:  ## constraint dummy point dist between start and end as zero otherwise infinity
                                    distances[from_counter][to_counter] = 0.0000000000000000000000001
                                else:
                                    distances[from_counter][to_counter] = 10000000000
                        else:
                            if to_counter == 0 or to_counter == len(locations)-2:  ## constraint dummy point dist between start and end as zero otherwise infinity
                                    distances[from_counter][to_counter] = 0.0000000000000000000000001
                            else:
                                distances[from_counter][to_counter] = 10000000000
            # ortools algorithm part
            manager = pywrapcp.RoutingIndexManager(len(data['location']), data['num_vehicles'], data['depot'])        
            rounting = pywrapcp.RoutingModel(manager)

            distance_matrix = distances   

            def distance_callback(from_index, to_index):
                from_node = manager.IndexToNode(from_index)
                to_node = manager.IndexToNode(to_index)
                return distance_matrix[from_node][to_node] 

            transit_callback_index = rounting.RegisterTransitCallback(distance_callback)

            rounting.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

            search_parameters = pywrapcp.DefaultRoutingSearchParameters()
            search_parameters.first_solution_strategy = (routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC)

            solution = rounting.SolveWithParameters(search_parameters) 

            if solution:
                index = rounting.Start(0)
                best_state = []
                route_distance = 0
                while not rounting.IsEnd(index):
                    best_state += [manager.IndexToNode(index)]
                    index = solution.Value(rounting.NextVar(index))

            # sorted to fixed start fixed stop sequence and filter dummy node out 
            n_state = len(best_state)

            ind_start = best_state.index(0)
            ind_end = best_state.index(n_state-2)
            ind_dummy = best_state.index(n_state-1)

            if ind_dummy > 0 and ind_dummy < n_state-1:
                if ind_start < ind_end:
                    ind_via_points = best_state[0:ind_start+1][::-1] + best_state[ind_end:][::-1]
                else:
                    ind_via_points = best_state[ind_start:] + best_state[0:ind_end+1]
            elif ind_dummy == 0:
                if ind_start == 1:
                    ind_via_points = best_state[1:]
                elif ind_start == n_state-1:
                    ind_via_points = best_state[1:][::-1]
            elif ind_dummy == n_state-1:
                if ind_start == n_state-2:
                    ind_via_points = best_state[0:n_state-1][::-1]
                elif ind_start == 0:
                    ind_via_points = best_state[0:n_state-1]

            # map index to via-points
            area_step = [locations[j] for j in ind_via_points]
            
            for i in range(len(area_step)-1):
                from_area = area_step[i]
                to_area = area_step[i+1]
                path_,dist = self.find_shortest_path(from_area,to_area)
                
                if i != len(area_step)-2:
                    path = path_[:-1]
                else:
                    path = path_

                for j in path:
                    self.area_path[r] += [j]
                    if j in self.cleanning_area[r]:
                        self.clean_path[r] += [1]
                        self.cleanning_area[r].pop(self.cleanning_area[r].index(j))
                    else:
                        self.clean_path[r] += [0]
            
            mid_edge_path = [list(set([e for e in self.area_graph._node[self.area_path[r][x]]['area_edges']]).intersection(set([e for e in self.area_graph._node[self.area_path[r][x+1]]['area_edges']]))) for x in range(len(self.area_path[r])-1)]
            
            filter_mid_edge_path = []
            for i in mid_edge_path:
                if len(i) > 1:
                    for j in i:
                        if not self.wall_graph[self.wall_edges_dict[j][0]][self.wall_edges_dict[j][1]]['isWall'] :
                            filter_mid_edge_path.append(j)
                else:
                    filter_mid_edge_path.append(i[0])

            self.mid_edge_path_idx[r] = [[self.mid_edges_graph._node[y]['edge_Number'] for y in self.mid_edges_graph._node].index(x)+1 for x in filter_mid_edge_path]

            print('finish generated area paths for robot : ' + str(r))
            print(self.area_path[r])
            print(self.clean_path[r])
            print('\n')

        return None
    
    def find_shortest_path(self, from_area, to_area):
        k = 1
        max_iter = self.max_iteration_area_search
        candidate = []
        current_graph_front = [from_area]
        current_graph_back = [to_area]
        connect_graph = []

        while True:
            while (k <= max_iter):
                temp_graph_front = []
                temp_graph_back = []
                        
                if type(current_graph_front[0]) == int: #first time boths are int
                    neighbor_list_front = [[n] for n in self.area_graph.neighbors(current_graph_front[0])]
                    temp_graph_front = np.concatenate((matlib.repmat(current_graph_front, len(neighbor_list_front),1), neighbor_list_front), axis = 1)

                    neighbor_list_back = [[n] for n in self.area_graph.neighbors(current_graph_back[0])]
                    temp_graph_back = np.concatenate((matlib.repmat(current_graph_back, len(neighbor_list_back),1), neighbor_list_back), axis = 1)
                else: 
                    graph_i_front = []
                    graph_i_back = []
                    num_front = len(current_graph_front[0])
                    num_back = len(current_graph_back[0])

                    for i in range(len(current_graph_front)):
                        neighbor_list_front = [[n] for n in self.area_graph.neighbors(current_graph_front[i][num_front-1])]
                        graph_i_front = np.concatenate((matlib.repmat(current_graph_front[i], len(neighbor_list_front), 1), neighbor_list_front), axis = 1)
                        
                        if i == 0:
                            temp_graph_front = graph_i_front
                        else:
                            temp_graph_front = np.concatenate((temp_graph_front, graph_i_front))
                        
                    for i in range(len(current_graph_back)):
                        neighbor_list_back = [[n] for n in self.area_graph.neighbors(current_graph_back[i][num_back-1])]
                        graph_i_back = np.concatenate((matlib.repmat(current_graph_back[i], len(neighbor_list_back), 1), neighbor_list_back), axis = 1)

                        if i == 0:
                            temp_graph_back = graph_i_back
                        else:
                            temp_graph_back = np.concatenate((temp_graph_back, graph_i_back))

                current_graph_front = temp_graph_front
                current_graph_back = temp_graph_back

                for i in range(len(current_graph_front)):
                    diff_list = np.subtract(current_graph_back, current_graph_front[i])
                    end_diff_list = np.transpose(diff_list)[len(diff_list[0])-1]
                    if current_graph_front[i][-1] == to_area:
                        candidate.append(current_graph_front[i])
                    if 0 in end_diff_list:
                        zero_index_list = [x for x in range(len(end_diff_list)) if end_diff_list[x] == 0]
                        [connect_graph.append(list(current_graph_front[i]) + list(current_graph_back[z][::-1][1:])) for z in zero_index_list]
                
                candidate += [path_list for path_list in connect_graph]
                
                k = k+1

            if len(candidate) == 0:
                max_iter = max_iter+1
            else:
                break
        
        all_dist = []
        for i in range(1, len(self.mid_edges_graph._node)+1):
            for j in self.mid_edges_graph[i]._atlas:
                all_dist.append(self.mid_edges_graph[i][j]['weight'])
                
        max_dist = max(all_dist)
        
        # find minimum path
        dist_path = []
        len_list = [len(x) for x in candidate]
        min_len = min(len_list)
        candi_i = candidate[len_list.index(min_len)]
        
        path_edge = [list(set([e for e in self.area_graph._node[candi_i[x]]['area_edges']]).intersection(set([e for e in self.area_graph._node[candi_i[x+1]]['area_edges']]))) for x in range(len(candi_i)-1)]
        filter_edge = []
        for i in path_edge:
            if len(i) > 1:
                for j in i:
                    if not self.wall_graph[self.wall_edges_dict[j][0]][self.wall_edges_dict[j][1]]['isWall'] :
                        filter_edge.append(j)
            else:
                filter_edge.append(i[0])

        node_idx = [[self.mid_edges_graph._node[y]['edge_Number'] for y in self.mid_edges_graph._node].index(x)+1 for x in filter_edge]
        dist_path= 0

        for idx in range(len(node_idx)-1):
            current_node = node_idx[idx]
            next_node = node_idx[idx+1]
            if next_node - current_node != 0:
                dist_path += self.mid_edges_graph[current_node][next_node]['weight']
                
        min_dist = dist_path
        min_path = list(candi_i)
        
        return min_path, min_dist
    
    def viaPointGenerator(self, vertices, current_pose, end_pose):
        [polygon, region_cm, region_radius] = self.coverage_points(vertices)
        region_cm = list(region_cm)
        region_radius = list(region_radius)

        via_points = []
        path_radius = []

        # create location datas for ortools
        data = {}
        data['location'] = [tuple(current_pose)]  # add start pose
        for i in region_cm:
            data['location'] += [tuple(i)] # add list of possibles pose
                
        dummy_point = (10000000000,10000000000)
        data['location'] += [tuple(end_pose), dummy_point] # add end pose and dummy node
        data['num_vehicles'] = 1
        data['depot'] = 1
        
        locations = data['location']
        
        # create distances matrix for ortools
        distances = {}
        for from_counter, from_node in enumerate(locations):
            distances[from_counter] = {}
            for to_counter, to_node in enumerate(locations):
                if from_counter == to_counter:
                    distances[from_counter][to_counter] = 0
                else:
                    if from_counter != len(locations)-1:
                        if to_counter != len(locations)-1:
                            distances[from_counter][to_counter] = (int(math.hypot((from_node[0] - to_node[0]),(from_node[1] - to_node[1]))))
                        else:
                            if from_counter == 0 or from_counter == len(locations)-2:  ## constraint dummy point dist between start and end as zero otherwise infinity
                                distances[from_counter][to_counter] = 0.0000000000000000000000001
                            else:
                                distances[from_counter][to_counter] = 10000000000
                    else:
                        if to_counter == 0 or to_counter == len(locations)-2:  ## constraint dummy point dist between start and end as zero otherwise infinity
                                distances[from_counter][to_counter] = 0.0000000000000000000000001
                        else:
                            distances[from_counter][to_counter] = 10000000000
        
        # ortools algorithm part
        manager = pywrapcp.RoutingIndexManager(len(data['location']), data['num_vehicles'], data['depot'])        
        rounting = pywrapcp.RoutingModel(manager)

        distance_matrix = distances   

        def distance_callback(from_index, to_index):
            from_node = manager.IndexToNode(from_index)
            to_node = manager.IndexToNode(to_index)
            return distance_matrix[from_node][to_node] 

        transit_callback_index = rounting.RegisterTransitCallback(distance_callback)

        rounting.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

        search_parameters = pywrapcp.DefaultRoutingSearchParameters()
        search_parameters.first_solution_strategy = (routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC)

        solution = rounting.SolveWithParameters(search_parameters) 

        if solution:
            index = rounting.Start(0)
            best_state = []
            route_distance = 0
            while not rounting.IsEnd(index):
                best_state += [manager.IndexToNode(index)]
                index = solution.Value(rounting.NextVar(index))

        # sorted to fixed start fixed stop sequence and filter dummy node out 
        n_state = len(best_state)

        ind_start = best_state.index(0)
        ind_end = best_state.index(n_state-2)
        ind_dummy = best_state.index(n_state-1)

        if ind_dummy > 0 and ind_dummy < n_state-1:
            if ind_start < ind_end:
                ind_via_points = best_state[0:ind_start+1][::-1] + best_state[ind_end:][::-1]
            else: 
                ind_via_points = best_state[ind_start:] + best_state[0:ind_end+1]
        elif ind_dummy == 0:
            if ind_start == 1:
                ind_via_points = best_state[1:]
            elif ind_start == n_state-1:
                ind_via_points = best_state[1:][::-1]
        elif ind_dummy == n_state-1:
            if ind_start == n_state-2:
                ind_via_points = best_state[0:n_state-1][::-1]
            elif ind_start == 0:
                ind_via_points = best_state[0:n_state-1]

        # map index to via-points
        for i in ind_via_points:
            if i != 0 and i != len(ind_via_points)-1:
                via_points.append(region_cm[i-1])
                path_radius.append(region_radius[i-1])
            
        return [via_points, path_radius]
    
    def coverage_points(self, vertices_i):
        radius = self.max_cleanning_radius
        cm = self.centroid(vertices_i)
        vertices = vertices_i[: -1]

        if any([i > radius**2 for i in [sum(np.subtract(ver,cm)**2) for ver in vertices]]):
            sub_poly = self.cutPolyMinor(vertices)
            s_1 = sub_poly[0]
            s_2 = sub_poly[1]
            [polygon_1, region_cm_1, region_radius_1] = self.coverage_points(s_1)
            [polygon_2, region_cm_2, region_radius_2] = self.coverage_points(s_2)
            cm_1 = self.centroid(sub_poly[0])
            cm_2 = self.centroid(sub_poly[1])
            
            if cm_1[0] < cm_2[0]:
                polygons = np.concatenate((polygon_1, polygon_2), axis = 1)
                region_cm = np.concatenate((region_cm_1, region_cm_2))
                region_radius = np.concatenate((region_radius_1, region_radius_2))
            else:
                polygons = np.concatenate((polygon_2, polygon_1), axis = 1)
                region_cm = np.concatenate((region_cm_2, region_cm_1))
                region_radius = np.concatenate((region_radius_2, region_radius_1))
                
            
        else:
            actual_radius_list = [np.sqrt(sum(np.subtract(ver,cm)**2)) for ver in vertices]
            idx = actual_radius_list.index(max(actual_radius_list))
            actual_radius = actual_radius_list[idx]
            polygons = [vertices]
            region_cm = [cm]
            region_radius = [actual_radius]

        return [polygons, region_cm, region_radius]

    def cutPolyMinor(self, vertices_i):
        vertices = np.transpose(vertices_i)
        cm = self.centroid(vertices_i)
        center = np.array([[cm[0]], [cm[1]]])
        dp = vertices - center
        xxyy = [sum(dp[0]**2), sum(dp[1]**2)]
        xy = sum(dp[0]*dp[1])
        I = [[xxyy[0], xy], [xy, xxyy[1]]]
        V = np.linalg.eigh(I)[1]
        theta_cut = -np.arctan2(V[0][0],V[1][0])
        rho_cut = np.cos(theta_cut)*center[0][0] + np.sin(theta_cut)*center[1][0]
        
        edge_temp = [ i for i in range(1, len(vertices[0])+1)]
        edge = np.transpose([edge_temp, edge_temp[1:] + edge_temp[:1]])

        sub_g = nx.Graph()
        sub_g_dict = {}

        # create nodes and edges
        for i in range(1, len(vertices[0])+1):
            sub_g.add_node(i, Coordinate = [vertices[0][i-1], vertices[1][i-1]] )
            sub_g.add_edge(edge[i-1][0], edge[i-1][1], label = i)
            sub_g_dict[i] = edge[i-1]

        num_ver = len(vertices[0])
        start_ver = 0
        next_idx = start_ver
        sub_poly = []
        for i in range(2):
            current_ver = next_idx
            loop_break = False
            coordinate = np.transpose(vertices)
            ver_list_minus = [coordinate[current_ver]]

            # search clock-wise until meets intersection
            while not loop_break:
                if current_ver != 0:
                    idx = current_ver - 1
                else:
                    idx = num_ver-1

                current = coordinate[current_ver]
                look_ahead = coordinate[idx]
                theta_look_ahead = -np.arctan2(look_ahead[0]-current[0], look_ahead[1]-current[1])
                rho_look_ahead = np.cos(theta_look_ahead)*current[0] + np.sin(theta_look_ahead)*current[1]

                A = [[np.cos(theta_look_ahead), np.sin(theta_look_ahead)], [np.cos(theta_cut), np.sin(theta_cut)]]
                b = [[rho_look_ahead], [rho_cut]]

                if (abs(theta_look_ahead-theta_cut) < 0.00000001) or (abs(abs(theta_cut-theta_look_ahead)-np.pi)<0.00000001):
                    ver_list_minus.append(look_ahead)
                    current_ver = idx
                else:
                    cut_intersect = np.matmul(np.linalg.inv(A),b)
                    if np.linalg.norm(cut_intersect-current) < 0.00000001:
                        loop_break = True
                    else:
                        cut_direction = -np.arctan2(cut_intersect[0]-current[0],cut_intersect[1]-current[1])
                        if np.cos(cut_direction-theta_look_ahead)<0:
                            ver_list_minus.append(look_ahead)
                            current_ver = idx
                        else:
                            distance_look_ahead = sum((look_ahead-current)**2)
                            distance_cut = sum((np.transpose(cut_intersect)-current)[0]**2)
                            if distance_cut <= distance_look_ahead:
                                ver_list_minus.append(np.transpose(cut_intersect)[0])
                                loop_break = True
                            else:
                                ver_list_minus.append(look_ahead)
                                current_ver = idx
            
            current_ver = next_idx
            loop_break = False
            coordinate = np.transpose(vertices)
            ver_list_plus = [coordinate[current_ver]]

            while not loop_break:
                if current_ver != num_ver-1:
                    idx = current_ver + 1
                else:
                    idx = 0

                current = coordinate[current_ver]
                look_ahead = coordinate[idx]
                theta_look_ahead = -np.arctan2(look_ahead[0]-current[0], look_ahead[1]-current[1])
                rho_look_ahead = np.cos(theta_look_ahead)*current[0] + np.sin(theta_look_ahead)*current[1]

                A = [[np.cos(theta_look_ahead), np.sin(theta_look_ahead)], [np.cos(theta_cut), np.sin(theta_cut)]]
                b = [[rho_look_ahead], [rho_cut]]

                if (abs(theta_look_ahead-theta_cut) < 0.00000001) or (abs(abs(theta_cut-theta_look_ahead)-np.pi)<0.00000001):
                    ver_list_plus.append(look_ahead)
                    current_ver = idx
                else:
                    cut_intersect = np.matmul(np.linalg.inv(A),b)
                    cut_direction = -np.arctan2(cut_intersect[0]-current[0],cut_intersect[1]-current[1]) 
                    if np.cos(cut_direction-theta_look_ahead)<0:
                        ver_list_plus.append(look_ahead)
                        current_ver = idx
                    else:
                        distance_look_ahead = sum((look_ahead-current)**2)
                        distance_cut = sum((np.transpose(cut_intersect)-current)[0]**2)
                        if (distance_cut <= distance_look_ahead) or (abs(distance_cut-distance_look_ahead) < 0.00000000001) and np.linalg.norm(cut_intersect-current) > 0.00000001:
                            ver_list_plus.append(np.transpose(cut_intersect)[0])
                            loop_break = True
                            if np.linalg.norm(distance_cut-distance_look_ahead) < 0.00000001:
                                if idx != num_ver-1:
                                    next_idx = idx+1
                                else:
                                    next_idx = 0
                            else:
                                next_idx = idx
                        else:
                            ver_list_plus.append(look_ahead)
                            current_ver = idx
        
            sub_poly.append(ver_list_minus[1:][::-1] + ver_list_plus)

        return sub_poly

    def calculate_time_kill(self, r):
        order = len(self.coef) - 1
        order_list = [ i for i in range(order,-1,-1)]
        time = 0

        if type(r) == list:
            for r_i in r:
                time_i = (self.energy_to_kill*100/(np.matmul((r_i ** order_list),self.coef)))[0]
                time += time_i
        else:
            time = (self.energy_to_kill*100/(np.matmul((r ** order_list),self.coef)))[0]
        
        return time

##############################################################################################################################
#####                                                                                                                    #####
#####                                            write to text file                                                      #####
#####                                                                                                                    #####
##############################################################################################################################

    def writeText(self):
        for r in range(len(self.robot_path)):
            file = open("robot "+str(r+1) + ".txt","w+")
            for i in range(len(self.robot_path[r])):
                if i==0:
                    continue
                t = 0
                if self.robot_duration[r][i] != 0:
                    t = 10
                qz, qw = yaw_to_qz_qw(self.robot_yaw[r][i]) 
                file.write(str(t)+','+str(self.robot_path[r][i][0])+','+str(self.robot_path[r][i][1])+','+str(qz)+','+str(qw)+'\n')
                #file.write(str(self.robot_duration[i])+','+str(self.robot_path[i][0])+','+str(self.robot_path[i][1])+','+str(self.robot_yaw[i])+','+str(self.robot_qw[i])+'\n')

##############################################################################################################################
#####                                                                                                                    #####
#####                                               visualization                                                        #####
#####                                                                                                                    #####
##############################################################################################################################

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

    def visualize_path(self):
        if self.VISUAL:
            self.fig = plt.figure()
            self.ax = plt.axes()
            self.line1, = self.ax.plot([],[], '--b', lw = 2)
            self.line2, = self.ax.plot([],[], 'r', lw = 2)

            for r in range(len(self.start_area)):
                if r == 0:
                    cl = 'k'
                elif r == 1:
                    cl = 'b'
                elif r == 2:
                    cl = 'r'

                plt.plot(self.K_data[self.K_idx==self.K_robot[r],0],self.K_data[self.K_idx==self.K_robot[r],1], 'o' + cl)
                plt.plot(self.K_centroids[self.K_robot[r],0],self.K_centroids[self.K_robot[r],1], 's' + cl, markersize=8)

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

            anim = animation.FuncAnimation(self.fig, self.animate, init_func = self.init_animate, frames = max([len(i) for i in self.robot_path]), interval = self.animate_times, blit = False)
            plt.axis('equal')
            plt.show()
    
    def init_animate(self):
        for r in range(len(self.robot_path)):
            self.line1.set_data([],[])
            self.line2.set_data([],[])
        
        return self.line1, self.line2

    def animate(self, i):
        l = 0.5
        for r in range(len(self.robot_path)):
            if r == 0:
                cl = 'k'
            elif r == 1:
                cl = 'b'
            elif r == 2:
                cl = 'r'
    
            if i <= len(self.robot_path[r])-1:
                p_i = self.robot_path[r][i]
                cir = plt.Circle((p_i[0], p_i[1]), self.robot_clean_radius[r][i], color = 'yellow')
                plt.plot(p_i[0], p_i[1], '*' + cl)
                self.ax.add_artist(cir)

                # plt.plot([p_i[0], p_i[0]+l*np.cos(self.robot_yaw[r][i])],[p_i[1], p_i[1]+l*np.sin(self.robot_yaw[r][i])], cl, lw = 2)
                
                if i != len(self.robot_path[r])-1:
                    p_i_p = self.robot_path[r][i+1]
                    plt.plot([p_i[0], p_i_p[0]], [p_i[1], p_i_p[1]], '--' + cl , lw = 2)
                    # plt.plot([p_i_p[0], p_i_p[0] + l*np.cos(self.robot_yaw[r][i+1])],[p_i_p[1], p_i_p[1] + l*np.sin(self.robot_yaw[r][i+1])], cl, lw = 2)
    
                if i != 0:
                    p_i_m = self.robot_path[r][i-1]
                    # plt.plot([[p_i_m][0], p_i[0]], [p_i_m[1], p_i[1]], '--k')

        return None
                
##############################################################################################################################
#####                                                                                                                    #####
#####                                              support function                                                      #####
#####                                                                                                                    #####
##############################################################################################################################

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

def yaw_to_qw(yaw):
    qw = np.cos(yaw/2)
    return qw

def yaw_to_qz_qw(yaw):
    qz = np.sin(yaw/2)
    qw = np.cos(yaw/2)
    return qz, qw 

def q_to_yaw(qz,qw):
    return math.atan2(2*qz*qw, 1-(2*qz*qz))

##############################################################################################################################
#####                                                                                                                    #####
##### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         test run Unit         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #####
#####                                                                                                                    #####
##############################################################################################################################

if __name__ == "__main__":

    #------------------------------------------------------------data from user interface--------------------------------------------------------------------------#
 
    # vertices (selected points)
    # x = [0, 2, -4, 2, -4,  2,  2, 0.5,  7, 14, 20.5, 22.5, 22.5, 16.5,   9,  4.5, 5, 8.5, 10.5, 7, 11.5, 14.5, 14, 18.5, 18.5, 16, 11.5, 11.5]
    # y = [0, 4,  6, 9, 11, 12, 14,  17, 18, 18, 14.5,  8.5,  2.5, -1.5, -1.5, 1.5, 6,   3,    5, 9,   10,  5.5,  2,    2,  7.5, 11,   14,   12]

    # x = [-75, -75, -75+110, -75+110, -75+110, -75+110, -75+110,  -75+110+45, -75+110+45, -75+110+45, -75+110+45, -75+110+50, -75+110+50, -75+110+50, -75+100+45+22, -75+100+45+22, -75+110+50+25, -75+110+50+25+45, -75+110+50+80, -75+110+192-40, -75+110+192-40, -75+110+192-40, -75+110+192, -75+110+192, -75+110+192, -75+110+192, -75+110+192+550, -75+110+192+550]
    # y = [10-40, 10, 10-40, 10, 10-40+100, 10-40+100+50, 10-40+100+50+562, 10-40+100+50, 10-40+100+50+45, 10-40+100+50+45+50, 10-40+100+50+562, 10-40, 10-40+60, 10-40+100, 10-40+100+50+45, 10-40+100+50+45+50, 10-40, 10-40-45, 10-40+60-50, 10-40-45, 10-40+60-50, 10-40+100, 10-40-45, 10-40, 10-40+100, 10-40+150, 10-40-45, 10-40]
    
    # x = np.array([-8, -4, 0, 2,  0,  -5,  5, 10,  3,   4,  6, 11, 17,  13, 21])
    # y = np.array([ 5, 0, 9, 5, -4, -10, -1,  2, -8, -15, -6, -8,  4, -15, -5])
    
    # x = [1.5, 1.5, 0.3, 0.3, -0.9, -0.9, -3.1, -3.1]
    # y = [0.9, -0.3, -0.3, 0.9, 0.3, -0.9, -0.9, 0.3]

    # x = [   -2.868809,   -3.178729,    2.134198,    2.458878,    1.263469,    2.665492,
    #         -2.662194,   -2.721227,   -2.824534,   -4.093734,   -5.23011,    -9.465694,
    #         -10.86772,   -11.88603,   -13.8341,    -13.71604,   -11.16288,   -12.29926,
    #         -10.86671,    -9.907436,   -9.376143,   -9.346628,   -7.487102,   -7.516618,
    #         -7.35428,    -7.295245,   -8.490656,   -8.667751,   -5.273382,   -5.450479,
    #         -4.181279,   -4.406027,  -14.1392,    -14.18191,   -15.1707,    -18.41749,
    #         -18.90451,   -14.78699,   -14.38835,   -13.57665,   -14.00463,    -8.898319,
    #         -7.378232,   -7.153488,  -12.59924,    -7.05362,    -4.515223,   -3.290298,
    #         -3.231264,   -3.54458,     4.882314,    5.3103,   ]
    # y =  [  1.006395,   -3.288222,   -3.657177,    0.5193763,   0.622683,    2.851161,
    #         3.279147,    2.70358,     1.640994,    1.109702,    1.19825,     1.508172,
    #         1.626236,    1.773818,    2.009948,    4.194151,    3.928505,   -3.893306,
    #         -3.947384,   -3.90311,    -3.932626,   -2.707702,   -2.825767,   -4.13924,
    #         -0.6858368,   0.6571531,   0.7899766,  -2.707702,    0.5243301,  -0.848177,
    #         -0.9219675,  -4.255771,   -2.338892,   -3.422133,   -2.639953,   -3.746814,
    #         -10.55346,   -10.87814,    -5.638999,   -5.72755,   -10.90765,   -11.23233,
    #         -10.39112,    -6.767299,   -5.719471,   -6.098875,   -5.877501,   -5.715162,
    #         -4.696852 , -10.0242 ,   -10.54074  ,  -4.371837 ]

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
    # p = np.array([x, y])*4/90

    # p = np.array([[-2, -1,  2,  2,  2, -1, -2, -1,  1],
    #               [-1, -1, -1,  1,  2,  2,  2,  1,  1]])

    # connectivity edges (drag and draw)
    # s =      [1, 1,  1, 2, 2,  2, 3, 3, 4, 4,  4, 5, 5, 6,  6,  6, 7,  7, 8,  9,  9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 17, 18, 19, 19, 19, 20, 20, 21, 21, 21, 22, 23, 24, 25, 26, 27]
    # t =      [2, 3, 16, 3, 4, 17, 4, 5, 5, 6, 20, 7, 8, 7, 20, 28, 8, 27, 9, 10, 27, 11, 27, 12, 26, 13, 25, 14, 24, 15, 23, 16, 18, 17, 18, 19, 20, 22, 23, 21, 28, 22, 26, 28, 23, 24, 25, 26, 27, 28]
    # isWall = [0, 1,  1, 0, 1,  1, 0, 1, 0, 0,  1, 0, 1, 1,  0,  1, 0,  1, 1,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  0,  1,  1,  1,  0,  0,  0,  0,  1,  1,  0,  1,  1,  1,  1,  0,  1]
    
    # s =      [1, 1, 2, 3, 3, 4, 5, 5,  6, 6, 7,  8, 8,  9,  9,  10, 10, 12, 12, 13, 13, 14, 15, 17, 17, 18, 19, 20, 20, 21, 22, 23, 23, 24, 24, 25, 27]
    # t =      [2, 3, 4, 4,12, 5, 6, 14, 7, 8, 11, 9, 26, 10, 15, 11, 16, 13, 17, 14, 19, 22, 16, 18, 19, 20, 21, 21, 23, 22, 25, 24, 27, 25, 28, 26, 28]
    # isWall = [1, 1, 1, 0, 1, 1, 1, 0,  1, 0, 1,  1, 1,  0,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,  0,  1,  1,  0,  1,  1,  0,  0,  1,  1,  1,  1,  1]

    # s = np.array([1, 1, 2, 2, 2, 3, 4, 4, 5, 5, 5,  6, 7,  7,  8,  8,  9,  9, 10, 10, 11, 12, 12, 13])
    # t = np.array([2, 3, 4 ,5 ,6 ,4, 7, 8, 6, 7, 9, 10, 8, 11, 12, 13, 10, 11, 12, 14, 12, 14, 15, 15])
    # isWall = np.array([1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1])

    # s = [1, 1, 2, 3, 3, 4, 5, 5, 6, 7]
    # t = [2, 4, 3, 4, 6, 5, 6, 8, 7, 8]
    # isWall = [0,0,0,0,0,0,0,0,0,0]

    # s       = [ 1 , 2,  3,  5,  1,  8,  7, 17, 16, 15, 33, 35, 36, 37, 38, 39, 39, 40, 41, 42, 43, 44, 46, 46,
    #             47, 48, 50, 51, 52,  2, 14, 18, 19, 13, 20, 21, 21, 22, 23, 24, 32, 31, 30, 25, 31, 12, 11, 29,
    #             26, 27,  5,  4,  8,  1, 10, 12, 14, 16, 33, 18, 40, 49, 19, 21, 28, 26, 29, 32]
    # t       = [ 2 , 3,  4,  1,  9,  7,  6,  7, 15, 33, 35, 36, 37, 38, 39, 34, 40, 41, 42, 43 ,44, 46, 45, 47,
    #             48, 50, 51, 52,  3, 49, 18, 19, 13, 14, 12, 20, 22, 28, 24, 32, 31, 30, 25, 23, 10, 11, 29, 26,
    #             27, 28,  4,  6,  9, 10, 11, 13, 15, 17, 34, 45, 45, 48, 20, 24, 23, 25, 30, 47]
    # isWall  = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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

    # selected_areas (consecutive_edges only)
    # area_edges = np.array([[1,4,2],[4,5,7],[7,9,8],[9,10,14,12],[12,17,13],[17,18,21,19],[21,23,20],[23,49,25,22],[25,48,27,24],[27,47,29,26],[29,46,31,28],[31,39,36,33,30],[33,35,34,32],[34,6,1,3],[39,45,38],[38,42,40,37],[40,44,41],[41,16,15],[11,15,10],[44,43,49,50]])
    # area_edges = np.array([[2, 4, 3, 1],[5, 18, 20, 8, 6, 4],[8, 22, 31, 36, 13, 10, 7],[10, 12, 14, 16, 11, 9],[15, 23, 17, 14],[19, 25, 21, 18],[24, 26, 28, 27, 25],[29, 32, 34, 31, 30, 28],[33, 37, 35, 32]])
    # area_edges = np.array([[1, 3, 6, 2], [3, 4, 10, 7], [4, 5, 9], [7, 13, 8], [9, 12, 17, 11], [13, 14, 21, 15], [17, 19, 21, 18], [19, 20, 22], [15, 23, 24, 16]])
    # area_edges = np.array([s[1, 2, 4, 3], [4, 6, 7, 5], [7, 8, 10, 9]])
    
    # area_edges = np.array( [[1, 2, 3, 51, 4], [6, 53, 5, 4, 51, 52, 7],
    #                         [9, 57, 34, 56, 46, 55, 54, 5, 53, 6, 8, 58],
    #                         [10, 59, 16, 17, 61, 60, 31, 57], [13, 14, 15, 16, 59, 11, 12],
    #                         [18, 19, 20, 21, 22, 23, 61],
    #                         [60, 23, 24, 68, 40, 64, 36, 63, 32], [33, 63, 35, 56],
    #                         [37, 64, 39, 65, 38], [50, 65, 44, 66, 49],
    #                         [66, 43, 67, 48], [47, 67, 42, 45, 55],
    #                         [45, 41, 68, 25, 62, 30, 1, 54], [30, 62, 26, 27, 28, 29, 2]])

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

    # connectivity areas (drag and draw)
    # s = [1,  1, 2, 3, 4,  4, 5, 6, 7, 8,  8,  9, 10, 11, 12, 12, 13, 15, 16, 17, 17, 18]
    # t = [2, 14, 3, 4, 5, 19, 6, 7, 8, 9, 20, 10, 11, 12, 13, 15, 14, 16, 17, 18, 20, 19]

    # s = [1, 2, 2, 3, 3, 4, 6, 7, 8]
    # t = [2, 3, 6, 4, 8, 5, 7, 8, 9]

    # s = [1, 2, 2, 3, 4, 5, 6, 6, 7]
    # t = [2, 3, 4, 5, 6, 7, 7, 9, 8]

    # s = [1, 2]
    # t = [2, 3]

    # s = [ 1,  2,  3,  3,  3,  3,  4,  4,  4, 13,  7,  7,  9, 10, 11,  7]
    # t = [ 2,  3, 13, 12,  8,  4,  5,  7,  6, 14,  8,  9, 10, 11, 12, 13]

    s = [37, 36, 25, 21, 21, 24, 24, 23, 22, 26, 28, 25, 27, 30, 27, 30, 30, 30, 30, 19,  1,  2,  2,  3,
    3,  3,  3,  3,  4,  4,  4, 13,  7,  7,  9, 10, 11,  7, 18, 17, 16, 15]
    t = [38, 37, 36, 22, 25, 25, 26, 24, 23, 28, 29, 27, 28, 31, 30, 32, 34, 33, 35, 20,  2, 21,  3, 13,
    12,  8,  4, 15,  5,  7,  6, 14,  8,  9, 10, 11, 12, 13, 19, 19, 20, 20]

    con_a = np.array([s, t])

    #----------------------------------Construction Graph---------------------------------------------------------------#
    CPP = coveragePathPlanning(p, con_e, isWall, area_edges, con_a, True)
    CPP.initialize_graphs()
    # CPP.visual_area()

    CPP.generate_robotWorkFlowPath([[24,-4],[0,0]], [[0,0],[24,-4]], [37,38,19,20,15,28,29,32,6])
    
    print("")