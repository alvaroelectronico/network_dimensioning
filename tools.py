import pandas as pd
import os
import time
from hard_coded_data import *

def read_data(folder_path):

# DATA_PATH = "..\..\..\datos_entrada\csv\casos_daniele"
# case_path = "0010km2_0"
# folder_path = os.path.join(DATA_PATH, case_path)

    start_time = time.time()
    df_initial_capacity = pd.read_csv(os.path.join(folder_path, "capacityi.csv"))
    df_potential_capacity = pd.read_csv(os.path.join(folder_path, "capacityp.csv"))
    df_existing_sites = pd.read_csv(os.path.join(folder_path, "existing_sites.csv"))
    df_potential_sites = pd.read_csv(os.path.join(folder_path, "potential_sites.csv"))
    df_demand = pd.read_csv(os.path.join(folder_path, "traffic_demand.csv"))
    df_coverage = pd.read_csv(os.path.join(folder_path, "coverage.csv"))
    print("csv files read: {}".format(time.time() - start_time))

    start_time = time.time()
    start_time = time.time()
    existing_sites = list(df_existing_sites.site_id.unique())
    potential_sites = list(df_potential_sites.site_id.unique())
    sites = existing_sites + potential_sites
    lots = df_demand.lot_id.unique()
    nodes = df_demand.node.unique()
    cells = df_coverage.cell.unique()

    existing_node_in_site = df_initial_capacity[df_initial_capacity.capacity > 0].drop_duplicates(['site_id', 'node']).set_index(['site_id', 'node']).index.to_list()
    potential_node_in_site = df_initial_capacity[df_initial_capacity.capacity == 0].drop_duplicates(['site_id', 'node']).set_index(['site_id', 'node']).index.to_list()
    existing_cell_in_site_node = df_initial_capacity[df_initial_capacity.capacity > 0].set_index(['site_id', 'node', 'cell']).index.to_list()
    potential_cell_in_site_node = df_initial_capacity[df_initial_capacity.capacity == 0].set_index(['site_id', 'node', 'cell']).index.to_list()
    print("existing and potential node-site and node-site-cell created: {}".format(time.time() - start_time))

    start_time = time.time()
    df_initial_capacity.set_index(['site_id', 'node', 'cell'], inplace=True)
    initial_capacity = df_initial_capacity.to_dict(orient="dict")['capacity']
    df_potential_capacity.set_index(['site_id', 'node', 'cell'], inplace=True)
    max_capacity = df_potential_capacity.to_dict(orient="dict")['capacity']
    df_demand.set_index(['lot_id', 'node'], inplace=True)
    demand = df_demand.to_dict(orient='dict')['demand']

    df_coverage.columns
    df_coverage2 = df_coverage.set_index(['site_id', 'node', 'cell', 'lot_id'])
    coverage = list(set(df_coverage2.index))
    print("coverage read: {}".format(time.time() - start_time))

    existing_node_in_site = list(set([(i[0], i[1]) for i in initial_capacity.keys() if initial_capacity[i] > 0]))

    potential_node_in_site = [(s, n) for s in sites for n in nodes
                              if not (s, n) in existing_node_in_site]

    existing_cell_in_site_node = [(i[0], i[1], i[2])
                                  for i in initial_capacity.keys() if initial_capacity[i] > 0]

    potential_cell_in_site_node = [i for i in initial_capacity.keys() if not i in existing_cell_in_site_node]

    start_time = time.time()
    df_coverage['site_cell'] = list(zip(df_coverage['site_id'], df_coverage['cell']))
    df_coverage3 = df_coverage.groupby(by=['lot_id', 'node'])['site_cell'].apply(lambda x: x.values.tolist())
    site_cells_lighting_lot_node = df_coverage3.to_dict()
    print("site_cells_lighting_lot_node read: {}".format(time.time() - start_time))

    start_time = time.time()
    #df_coverage['site_node_cell'] = list(zip(df_coverage['site_id'], df_coverage['node'], df_coverage['cell']))
    df_coverage4 = df_coverage.groupby(by=['site_id', 'node', 'cell'])['lot_id'].apply(lambda x: x.values.tolist())
    lots_covered_by_site_node_cell = df_coverage4.to_dict()
    print("lots_covered_by_cells read: {}".format(time.time() - start_time))

    return lots, sites, nodes, cells, existing_sites, potential_sites, initial_capacity, max_capacity, demand, \
           coverage, existing_node_in_site, potential_node_in_site, existing_cell_in_site_node, potential_cell_in_site_node, \
           site_cells_lighting_lot_node, lots_covered_by_site_node_cell


def check_solution(self):
    if len(self.solution.keys()) == 0:
        self.get_solution()

    # Demand is met
    # demand_unmet = [(l, n) for l in self.lots for n in self.nodes
    #                 if self.demand[l, n] > sum(self.solution['traffic_of_cell'][s, n, c, l]
    #                                           for s in self.sites for c in self.cells
    #                                           if (s, n, c, l) in self.coverage)
    #                ]
    # if len(demand_unmet) > 0:
    #     for (l, n) in demand_unmet:
    #         self.solution_check += "ERROR. Demand not met. Lot {}, node {}\n".format(l, n)
    # else:
    #      self.solution_check += "OK. Demand constraint met"

    # Capacity not violated
    capacity_violated = [(s, n, c) for s in self.sites for n in self.nodes for c in self.cells
                         if self.solution["final_capacity"][s, n, c] < sum(self.solution['traffic_of_cell'][s, n, c, l]
                                                                           for l in self.lots if
                                                                           (s, n, c, l) in self.coverage)]
    if len(capacity_violated) > 0:
        for (s, n, c) in capacity_violated:
            self.solution_check += "ERROR. Capacity violated. Site {}, node {}, cell {}\n".format(s, n, c)
    else:
        self.solution_check += "OK. Capacity not violated"
    print(self.solution_check)

def export_sol_to_csv(solution, ouput_path, file_name):
    # Traffic
    df = pd.Series(solution["traffic"]).reset_index()
    df.columns = ["site", "node", "cell", "lot", "traffic" ]
    df.to_csv(os.path.join(ouput_path, "traffic_{}.csv".format(file_name)))
    # New sites
    if len(solution["new_sites"]) > 0:
        df = pd.Series(solution["new_sites"])
        df.columns = ["site"]
        df.to_csv(os.path.join(ouput_path, "new_sites_{}.csv".format(file_name)))
    # New nodes
    if len(solution["new_nodes"]) > 0:
        df = pd.DataFrame.from_dict(solution["new_nodes"], orient="columns")
        df.columns = ["site", "node"]
        df.to_csv(os.path.join(ouput_path,"new_nodes_{}.csv".format(file_name)))
    # New cells
    if len(solution["new_cells"]) > 0:
        df = pd.DataFrame.from_dict(solution["new_cells"], orient="columns")
        df.columns = ["site", "node", "cell"]
        df.to_csv(os.path.join(ouput_path,"new_cells_{}.csv".format(file_name)))
    # Upgraded cells
    if len(solution["upgraded_cells"]) > 0:
        df = pd.DataFrame.from_dict(solution["upgraded_cells"], orient="columns")
        df.columns = ["site", "node", "cell"]
        df.to_csv(os.path.join(ouput_path, "upgraded_cells_{}.csv".format(file_name)))
    # pd.DataFrame(performance_data, index=[0]).to_csv(os.path.join(ouput_path, "general_{}.csv".format(file_name)))

def recalculate_obj_func(solution):
    cost = len(solution['new_sites'])*pCAPEX_NEW_SITE
    cost += len(solution['new_nodes'])*pCAPEX_NEW_NODE
    cost += len(solution['new_cells'])*pCAPEX_NEW_CELL
    cost += len(solution['upgraded_cells'])*pCAPEX_UPGRADE_CEll
    cost += len(solution['new_sites'])*pOPEX_SITE
    return cost


def read_output_daniele(file_path):
 df = pd.read_csv(file_path)
 return df

def solution_check(df):
    traffic = df.set_index(['Site', 'node', 'cell', 'lot'])["demand"].to_dict()
    # capacity = df.set_index(['Site', 'node', 'cell'])['capacity'].to_dict()

    df_max_cap_violated = df[df.capacity > df.max_cap]

    df_traffic_cell = df.groupby(by=["Site", "node", "cell"])["demand"].sum().reset_index()
    df_traffic_cell = pd.merge(df_traffic_cell, df[["Site", "node", "cell", "capacity"]], how="left")
    df_cap_violated = df_traffic_cell[df.demand > df.capacity]

    df_demand_satified = df.groupby(by=["node", "lot"])["demand"].sum().reset_index()

    pass
