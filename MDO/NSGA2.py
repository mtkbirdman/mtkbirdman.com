from CST import CST # 自作

import sys
import time
import numpy as np
import pandas as pd

from xfoil import XFoil
from xfoil.model import Airfoil
from ctypes import cdll

from pymoo.util.misc import stack
from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.sampling.lhs import LHS
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PolynomialMutation
#from pymoo.termination import get_termination # pymoo==0.6.0
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter

# Create an instance of the XFoil class for airfoil analysis
xf = XFoil()

# Define a custom optimization problem
class MyProblem(Problem):
    def __init__(self, file_name='airfoil_shape.csv', n_wu=4, n_wl=4, x_diff=0.8):
        self.n_wu = n_wu  # Number of control points for the upper surface
        self.n_wl = n_wl  # Number of control points for the lower surface
        
        # Create an instance of the CST class for airfoil parameterization
        airfoil = CST()
        airfoil.fit_CST(file_name=file_name, n_wu=self.n_wu, n_wl=self.n_wl)
        
        # Retrieve airfoil properties
        self.thickness, self.x_tmax = airfoil.get_thickness()
        self.ler = airfoil.get_ler()
        self.TE_angle = airfoil.get_TE_angle()
        self.section_area = airfoil.get_section_area()
        self.w_LE = np.abs(np.abs(airfoil.wu[0,0])-np.abs(airfoil.wl[0,0]))
        
        # Initial guess for the design variables
        x0 = np.concatenate([airfoil.wu, airfoil.wl], axis=1).reshape(-1)
        xu = x0*(1 + np.sign(x0)*x_diff)  # Upper bounds for design variables
        xl = x0*(1 - np.sign(x0)*x_diff)  # Lower bounds for design variables

        super().__init__(n_var=self.n_wu + self.n_wl,
                         n_obj=2,
                         n_constr=6,
                         xl=xl,
                         xu=xu)
        
        # Additional properties for airfoil analysis
        self.M = 0.0  # Mach number
        self.n_cirt = 9  # Critical amplification ratio
        self.xtr = (1, 1)  # Transition points on the upper and lower surfaces
        self.max_iter = 100  # Maximum number of iterations for XFoil
        self.Re = 1000000  # Reynolds number
        self.cl = 0.5  # Target lift coefficient

    def _evaluate(self, X, out, *args, **kwargs):
        wu = X[:, 0:self.n_wu].tolist()  # Upper surface control points
        wl = X[:, self.n_wu:].tolist()  # Lower surface control points
        dz = np.zeros(X.shape[0]).reshape(-1, 1).tolist()  # Zero thickness distribution
        Node = 201  # Number of nodes for airfoil discretization
        
        # Create airfoil shapes using CST
        airfoil = CST()
        airfoil.create_airfoil(wu=wu, wl=wl, dz=dz, Node=Node)
        
        # Retrieve updated airfoil properties
        thickness, x_tmax = airfoil.get_thickness()
        ler = airfoil.get_ler()
        TE_angle = airfoil.get_TE_angle()
        section_area = airfoil.get_section_area()
        w_LE = np.abs(np.abs(airfoil.wu[:,0])-np.abs(airfoil.wl[:,0]))
        
        f1_list = []  # List to store drag coefficients
        f2_list = []  # List to store pitching moments
        g1_list = []  # List to store lift coefficients (for constraints)
        
        for i in range(X.shape[0]):
            # Calculate aerodynamic properties for each airfoil shape
            f1, f2, g1 = self._calculate(airfoil=airfoil, idx=i)
            f1_list.append(f1)
            f2_list.append(f2)
            g1_list.append(g1)
        
        f1 = np.array(f1_list)
        f2 = np.array(f2_list)

        g_diff = 0.5
        g1 = np.array(g1_list)  # cl0 > 0
        g2 = (self.thickness - thickness) / self.thickness
        g3 = (w_LE-self.w_LE/g_diff) / self.w_LE
        g4 = (self.section_area*g_diff - section_area) / self.section_area
        g5 = (self.TE_angle*g_diff - TE_angle) / self.TE_angle
        g6 = (self.ler*g_diff - ler) / self.ler
        #g7 = (0.00 - x_tmax) / self.x_tmax
        
        # Output objectives and constraints
        out["F"] = np.column_stack([f1, f2])
        out["G"] = np.column_stack([g1, g2, g3, g4, g5, g6]) #, g7])

    def _calculate(self, airfoil, idx):
        # Assign airfoil coordinates to XFoil
        xf.airfoil = Airfoil(airfoil.x[idx], airfoil.y[idx])
        
        # Define analysis parameters
        xf.M = self.M  # Mach number
        xf.n_crit = self.n_cirt  # Critical amplification ratio
        xf.xtr = self.xtr  # Transition points
        xf.max_iter = self.max_iter  # Maximum number of iterations
        xf.Re = self.Re  # Reynolds number
        
        # Perform analysis with XFoil
        xf.repanel()  # Re-panel airfoil
        a, cd, cm, cp = xf.cl(self.cl)  # Analyze at fixed lift coefficient
        cl0, cd0, cm0, cp0 = xf.a(0)  # Analyze at zero angle of attack
        xf.reset_bls()  # Reset boundary layers
        
        return cd, -cm, -cl0  # Return drag, pitching moment, and lift coefficient

if __name__ == '__main__':
    # Mesure execution time
    start_time = time.time()

    # Define the optimization problem
    problem = MyProblem(file_name='./MDO/NACA4412.csv', n_wu=5, n_wl=5, x_diff=1.5)
    
    # Initialize the NSGA-II algorithm
    algorithm = NSGA2(
        pop_size=100,
        #n_offsprings=10,
        sampling=LHS(),  # Latin Hypercube Sampling
        crossover=SBX(prob=0.9, eta=15),  # Simulated Binary Crossover
        mutation=PolynomialMutation(eta=20),  # Polynomial Mutation
        eliminate_duplicates=True
    )

    # Termination criterion (optional)
    #termination = get_termination("n_gen", 40) # pymoo==0.6.0
    
    # Execute the optimization
    res = minimize(
        problem,
        algorithm,
        #termination, # pymoo==0.6.0
        ("n_gen", 100),  # Number of generations
        seed=1,  # Random seed for reproducibility
        save_history=True,
        verbose=True
    )

    # Save optimization results
    pd.DataFrame(res.X).to_csv('./MDO/res_X.csv', index=False)
    pd.DataFrame(res.F, columns=['cd', 'cm']).to_csv('./MDO/res_F.csv', index=False)
    pd.DataFrame(res.G, columns=['cl0', 't/c', 'w_LE', 'section_area', 'TE_angle', 'ler']).to_csv('./MDO/res_G.csv', index=False)
    #pd.DataFrame(res.G, columns=['cl0', 't/c', 'w_LE', 'section_area', 'TE_angle', 'ler', 'x_tmax', ]).to_csv('./MDO/res_G.csv', index=False)
    
    # Plot the objective space
    plot = Scatter(title="Objective Space")
    plot.add(res.F)
    plot.show()

    end_time = time.time()
    print(f"Run time: {end_time - start_time} s")
