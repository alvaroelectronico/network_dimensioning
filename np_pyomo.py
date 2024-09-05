import pandas as pd
import os
from pyomo.environ import *
import time
from tools import *
from hard_coded_data import *
import pickle

# To debug this py file (to be commented or deleted):
DATA_PATH = "..\..\..\datos_entrada\csv\casos_20220401"
case_path = "00001km2_0"
# case_path = "010_12scope_2020_13_1805"
# case_path = "001_1"
folder_path = os.path.join(DATA_PATH, case_path)
print("Reading data")
start_time = time.time()

class NP_pyomo:
    def __init__(self, name, input_folder):
        self.name = name
        self.input_folder = input_folder
        self.lots = list()
        self.sites = list()
        self.nodes = list()
        self.cells = list()
        self.existing_sites = list()
        self.potential_sites = list()
        self.initial_capacity = dict()
        self.max_capacity = dict()
        self.demand = dict()
        self.coverage = list()
        self.existing_node_in_site = list()
        self.potential_node_in_site = list()
        self.existing_cell_in_site_node = list()
        self.potential_cell_in_site_node = list()
        self.site_cells_lighting_lot_node = list()
        self.lots_covered_by_site_node_cell = list()
        self.solver_params = dict()
        self.read_data()
        self.errors = dict()
        self.solution = dict()
        self.data = dict()
        self.model_instance = ConcreteModel()

    def read_data(self):
        self.lots, self.sites, self.nodes, self.cells, self.existing_sites, self.potential_sites, self.initial_capacity,\
        self.potential_capacity, self.demand, self.coverage, self.existing_node_in_site, self.potential_node_in_site,\
        self.existing_cell_in_site_node, self.potential_cell_in_site_note, self.site_cells_lighting_lot_node, \
        self.lots_covered_by_site_node_cell \
        = read_data(self.input_folder)
        print("Data read. Time: {}".format(time.time() - start_time))

    def prepare_data(self):
        # Imput data
        data = {None: {
            'sLots': {None: self.lots},
            'sSites': {None:  self.sites},
            'sNodes': {None:  self.nodes},
            'sCells': {None:  self.cells},
            'sCoverage': {None:  self.coverage},
            'sLotsNodes': {None: list( self.demand.keys())},
            'sExistingSites': {None:  self.existing_sites},
            'sPotentialSites': {None:  self.potential_sites},
            'sPotentialNodesInSites': {None:  self.potential_node_in_site},
            'sExistingNodesInSites': {None:  self.existing_node_in_site},
            'sExistingCellsInSitesNodes': {None:  self.existing_cell_in_site_node},
            'sPotentialCellsInSitesNodes': {None:  self.potential_cell_in_site_note},
            'pTrafficDemand':  self.demand,
            'pInitialCapacity':  self.initial_capacity,
            'pMaxCapacity':  self.potential_capacity
        }}
        return data

    def create_model(self):
        # Definint the abstract model
        model = AbstractModel()

        # Sets definition
        model.sSites = Set()
        model.sExistingSites = Set()
        model.sPotentialSites = Set()
        model.sNodes = Set()
        model.sLots = Set()
        model.sCells = Set()

        # Compound sets
        model.sCoverage = Set(dimen=4)
        model.sLotsNodes = Set(dimen=2)
        model.sExistingNodesInSites = Set(dimen=2)
        model.sPotentialNodesInSites = Set(dimen=2)
        model.sPotentialCellsInSitesNodes = Set(dimen=3)
        model.sExistingCellsInSitesNodes = Set(dimen=3)

        # Parameters
        model.pTrafficDemand = Param(model.sLots, model.sNodes, mutable=True)
        model.pInitialCapacity = Param(model.sSites, model.sNodes, model.sCells, mutable=True)
        model.pMaxCapacity = Param(model.sSites, model.sNodes, model.sCells, mutable=True)
        model.pCapexUpgradeCell = Param(mutable=True)
        model.pCapexNewCell = Param(mutable=True)
        model.pCapexNewNode = Param(mutable=True)
        model.pCapexNewSite = Param(mutable=True)
        model.pOpexSite = Param(mutable=True)
        model.pOpexNode = Param(mutable=True)

        # Location variables
        model.v01NewSite = Var(model.sPotentialSites, domain=Binary)
        model.v01NewNode = Var(model.sPotentialNodesInSites, domain=Binary)
        model.v01NewCell = Var(model.sPotentialCellsInSitesNodes, domain=Binary)
        model.v01UpgradeCell = Var(model.sExistingCellsInSitesNodes, domain=Binary)
        model.vFinalCapacity = Var(model.sSites, model.sNodes, model.sCells, domain=NonNegativeReals)

        # Assignment variables
        model.vTrafficOfCell = Var(model.sSites, model.sNodes, model.sCells, model.sLots, domain=NonNegativeReals)
        model.vTotalCost = Var(domain=NonNegativeReals)


        # Constraint rules
        def fcDemandFulfilment(model, lot, node):
            return sum(model.vTrafficOfCell[site, node, cell, lot]
                       for site in model.sSites for cell in model.sCells
                       if (site, node, cell, lot) in model.sCoverage) >= model.pTrafficDemand[lot, node]


        def fcMinCellCapacity(model, site, node, cell):
            return model.vFinalCapacity[site, node, cell] >= model.pInitialCapacity[site, node, cell]


        def fcMaxCellCapacity(model, site, node, cell):
            if (site, node, cell) in model.sExistingCellsInSitesNodes:
                return model.vFinalCapacity[site, node, cell] <= \
                       model.pInitialCapacity[site, node, cell] * (1 - model.v01UpgradeCell[site, node, cell]) \
                       + model.pMaxCapacity[site, node, cell] * model.v01UpgradeCell[site, node, cell]
            else:
                return model.vFinalCapacity[site, node, cell] <= \
                       model.pMaxCapacity[site, node, cell] * model.v01NewCell[site, node, cell]


        def fcNewCellIfNodeExists(model, site, node, cell):
            if (site, node) in model.sExistingNodesInSites:
                return Constraint.Skip
            else:
                return model.v01NewCell[site, node, cell] <= model.v01NewNode[site, node]


        def fcNewNodeIfSiteExists(model, site, node):
            if site in model.sExistingSites:
                return Constraint.Skip
            else:
                return model.v01NewNode[site, node] <= model.v01NewSite[site]


        def fcMaxCellTraffic(model, site, node, cell):
            return sum(model.vTrafficOfCell[site, node, cell, lot] for lot in model.sLots if
                       (site, node, cell, lot) in model.sCoverage) <= model.vFinalCapacity[site, node, cell]

        def fcEnoughCapacityGlobal(model, node):
            return sum(model.vFinalCapacity[site, node, cell] for site in model.sSites for node in model.sNodes
                       for cell in model.sCells) >= sum(model.pTrafficDemand[lot, node] for lot in model.sLots)

        def fcEnoughCapacityPerLot(model, node, lot):
            return sum(model.vFinalCapacity[site, node, cell] for site in model.sSites for node in model.sNodes
                       for cell in model.sCells if (site, node, cell, lot) in model.sCoverage) \
                   >= model.pTrafficDemand[lot, node]

        # Total cost
        def fvTotalCost(model):
            return model.vTotalCost == \
                   pCAPEX_NEW_SITE * sum(model.v01NewSite[site] for site in model.sPotentialSites) \
                   + pCAPEX_NEW_NODE * sum(model.v01NewNode[site, node] for (site, node) in model.sPotentialNodesInSites) \
                   + pCAPEX_NEW_CELL * sum(
                model.v01NewCell[site, node, cell] for (site, node, cell) in model.sPotentialCellsInSitesNodes) \
                   + pCAPEX_UPGRADE_CEll * sum(
                model.v01UpgradeCell[site, node, cell] for (site, node, cell) in model.sExistingCellsInSitesNodes) \
                   + pOPEX_SITE * sum(model.v01NewSite[site] for site in model.sPotentialSites) \
                   + pOPEX_NODE + sum(model.v01NewNode[site, node] for (site, node) in model.sPotentialNodesInSites)


        # Objective function
        def obj_expression(model):
            return model.vTotalCost


        # def activate_constratins(model):
        # Activating constraints

        # model.cDemandFulfilment = Constraint(model.sLotsNodes, rule=fcDemandFulfilment)
        model.cMinCellCapacity = Constraint(model.sExistingCellsInSitesNodes, rule=fcMinCellCapacity)
        model.cMaxCellCapacity = Constraint(model.sSites, model.sNodes, model.sCells, rule=fcMaxCellCapacity)
        model.cMaxCellTraffic = Constraint(model.sSites, model.sNodes, model.sCells, rule=fcMaxCellTraffic)
        model.cNewCellIfNodeExists = Constraint(model.sPotentialCellsInSitesNodes, rule=fcNewCellIfNodeExists)
        model.cNewNodeIfSiteExists = Constraint(model.sPotentialNodesInSites, rule=fcNewNodeIfSiteExists)
        model.cvTotalCost = Constraint(rule=fvTotalCost)
        model.cEnoughCapacityGlobal = Constraint(model.sNodes, rule=fcEnoughCapacityGlobal)
        model.cEnoughCapacityPerLot = Constraint(model.sNodes, model.sLots, rule=fcEnoughCapacityPerLot)

        # Objective function
        model.obj_func = Objective(rule=obj_expression)

        return model

    def get_solution(self):
        self.solution['new_sites'] = [s for s in self.potential_sites if self.model_instance.v01NewSite[s].value == 1]
        self.solution['new_nodes'] = [(s, n) for (s, n) in self.potential_node_in_site if self.model_instance.v01NewNode[s, n].value == 1]
        self.solution['new_cells'] = [(s, n, c) for (s, n, c) in self.potential_cell_in_site_node if
                                      self.model_instance.v01NewCell[s, n, c].value == 1]
        self.solution['upgraded_cells'] = [(s, n, c) for (s, n, c) in self.existing_cell_in_site_node if
                                           self.model_instance.v01UpgradeCell[s, n, c].value == 1]
        self.solution['traffic'] = {i: self.model_instance.vTrafficOfCell[i].value for i in self.coverage}
        self.solution['final_capacity'] = {(s, n, c): self.model_instance.vFinalCapacity[s, n, c].value for s in self.sites for n in
                                           self.nodes for c in self.cells}

    def set_solver_params(self, opt):
        for param in self.solver_params.keys():
            opt.options[param] = self.solver_params[param]

    def solve_model(self):
        print("Preparing data")
        data = self.prepare_data()
        print("Creating abstract model")
        model = self.create_model()
        print("Creating model instance")
        self.model_instance = model.create_instance(data)
        opt = SolverFactory('gurobi')
        self.set_solver_params(opt)
        print("Solving model instance")
        results = opt.solve(self.model_instance)
        print("{} solved\n".format(self.name))

    def print_info(self):
        pass
        # Creating model instance
        # print("Building model")
        # start_time = time.time()
        # instance = model.create_instance(input_data)
        # print("Model build. Time: {}".format(time.time() - start_time))
        #
        # # Setting the solver
        # opt = SolverFactory('gurobi')
        # print("Solving")
        # results = opt.solve(instance, tee=True)
        # print("Solving. Time: {}".format(time.time() - start_time))
        #
        # print("Total cost: {}".format(instance.vTotalCost.value))
        #
        # print()
        # print("New sites")
        # for site in [s for s in instance.sPotentialSites if instance.v01NewSite[s].value == 1]:
        #     print("{}".format(site))
        #
        # print()
        # print("New nodes")
        # for (site, node) in [(s, n) for (s, n) in instance.sPotentialNodesInSites if instance.v01NewNode[s, n].value == 1]:
        #     print("{}, {}".format(site, node))
        #
        # print()
        # print("New cells")
        # for (site, node, cell) in [(s, n, c) for (s, n, c) in instance.sPotentialCellsInSitesNodes if instance.v01NewCell[s, n, c].value == 1]:
        #     print("{}, {}, {}".format(site, node, cell))
        #
        #
        # print()
        # print("Upgraded cells")
        # for (site, node, cell) in [(s, n, c) for (s, n, c) in instance.sExistingCellsInSitesNodes if instance.v01UpgradeCell[s, n, c].value == 1]:
        #     print("{}, {}, {}".format(site, node, cell))
        #
        # print("Total cost: {}".format(instance.vTotalCost.value))

def main():
    # DATA_PATH = "C:\\Users\Alvaro\Dropbox\\academico\\network_planning\instances"
    DATA_PATH = "D:\Dropbox\\academico\\network_planning\instances" # Dropbox Lenovo

    cases_paths = {c: os.path.join(DATA_PATH, c) for c in os.listdir(DATA_PATH) if
                   os.path.isdir(os.path.join(DATA_PATH, c))}

    df_performance = pd.DataFrame(columns=["case", "obj_func", "gap", "run_time"])

    ordered_cases_paths_keys = list(cases_paths.keys())
    ordered_cases_paths_keys.sort()

    for case in ordered_cases_paths_keys[:2]:
        instance = NP_pyomo(case, cases_paths[case])
        instance.solver_params = dict(TIME_LIMIT=6000, MIPGap=0.00)
        instance.solve_model()
        instance.get_solution()
        # instance.get_df_performance_data()
        # df_performance = df_performance.append(instance.performance_data, ignore_index=True)
        # df_performance.to_csv(os.path.join(DATA_PATH, "results_network_planning.csv"))

        recalculated_cost = recalculate_obj_func(instance.solution)
        export_sol_to_csv(instance.solution, cases_paths[case], instance.name)

        with open(os.path.join(DATA_PATH, cases_paths[case], "solution.pickle"), 'wb') as handle:
            pickle.dump(instance.solution, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print("case {} solved".format(case))

    # df_performance.to_csv("network_planning.csv")

if __name__ == "__main__":
    main()
    print("run completed")