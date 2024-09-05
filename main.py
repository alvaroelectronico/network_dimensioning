import pandas as pd
from hard_coded_data import *
import gurobipy as gp
import os
from gurobipy import GRB
from tools import read_data
import time
import pickle
from np_gurobipy import NP_gurobipy
from tools import *

def main():
    # DATA_PATH = "C:\\Users\Alvaro\Dropbox\\academico\\network_planning\instances"
    DATA_PATH = "D:\Dropbox\\academico\\network_planning\instances" # Dropbox Lenovo

    cases_paths = {c: os.path.join(DATA_PATH, c) for c in os.listdir(DATA_PATH) if os.path.isdir(os.path.join(DATA_PATH, c))}

    df_performance = pd.DataFrame(columns=["case", "obj_func", "gap", "run_time"])

    ordered_cases_paths_keys = list(cases_paths.keys())
    ordered_cases_paths_keys.sort()

    for case in ordered_cases_paths_keys[:2]:
        instance = NP_gurobipy(case, cases_paths[case])
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
