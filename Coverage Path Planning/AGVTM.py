import numpy as np
import networkx as nx
import random
import matplotlib.pyplot as plt

class agv:
    
    def __init__(self, agv_specs):
        self.max_step = agv_specs[0]
        self.max_gap = agv_specs[1]
        self.max_incline = agv_specs[2]
        self.max_payload = agv_specs[3]
        self.velocity = agv_specs[4]




class constructionGraph(agv):
    graph = nx.Graph()
    sigma = 0.5
    max_distance = 1000
    
    def __init__(self, num_node, agv_specs):
        super().__init__(agv_specs)
        n = num_node
        self.num_node = n
        
        ## create nodes
        self.node_labels = {}
        for node in range(num_node):
            self.node_labels[node+1] = str(node+1)
            self.graph.add_node(node+1)
            
        ## create edges
        # random edge
        edge_dicts = {}
        for edge in range(int(num_node*2.5)):
            edge_dicts[edge+1] = [random.randint(1,num_node), random.randint(1,num_node)]
            
            
        # add edge
        for edge in edge_dicts:
            #random edge properties
            step = random.uniform(0,1+self.sigma)*self.max_step
            gap = random.uniform(0,1+self.sigma)*self.max_gap
            incline = random.uniform(0,1+self.sigma)*self.max_incline
            distance = random.uniform(0,1+self.sigma)*self.max_distance
            
            from_node = edge_dicts[edge][0]
            to_node = edge_dicts[edge][1]
            
            self.graph.add_edge(from_node, to_node, step= step, gap= gap, incline = incline, distance = distance)
    
    def visualize_graph(self):
        nx.draw(self.graph, labels = self.node_labels)
        plt.show()


    
class AGVTP(constructionGraph):
    
    def __init__(self, num_node = 20, agv_specs = [15, 20, 5, 250, 200], num_loads = 15):
        super().__init__(num_node, agv_specs)
       
        self.load_list = []    
        #random payload 
        for i in range(num_loads):
            load = random.uniform(0,1+self.sigma)*self.max_payload
            self.load_list.append(load)
            
        self.job_list = []
        # random job_list
        for num_job in range(num_loads):
            self.job_list.append([random.randint(1,num_node), random.randint(1,num_node)])
        
        self.job_list = np.transpose(self.job_list).tolist()

    def generate_path(self, load_list, location_list, start_node, goal_node):
        
        ###### develop code here #############
        
        pass
        
if __name__ == "__main__":
    num_nodes = 40
    agv_specs = [15, 20, 5, 250, 200]
    num_loads = 15
    
    start_node = 1
    goal_node = 15
    
    agv_planner = AGVTP(num_nodes, agv_specs, num_loads)
    agv_planner.visualize_graph()
    agv_planner.generate_path(agv_planner.load_list, agv_planner.job_list, start_node, goal_node)