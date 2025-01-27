from system_builder import SystemBuilder
from matplotlib import pyplot as plt
import numpy as np
import parameters as pr
import affichage as aff

""" 
This is the main script to run the FDM.

There are 2 ways to use this script:
1. The user can change the parameters of the system_builder at init_params
to solve different problems.
2. The user can use the function traitement() from the affichage.py to
choose the parameters of the system_builder.

To change from one way to another, the user must comment the other way.

Also, the user can choose which solution to display by using the functions
load_2_canva_diff_solution(), load_2_canva_numerical_solution(),
load_2_canva_analytical_solution() and load_combined_canvas() from the system_builder.py.

Important notes:
- If the user want to implement the FDM for 
the final problem, make sure to set the 
"final_setup" parameter to True in the multiplier_4_increasing_resolution
This parameters is not used because it will used the image to find 
the resolution of the grid given Lx and Ly.
"""


init_params = {
    "Lx": 1,
    "Ly": 1,
    "bc_side": 0,               # 0: "dirichlet" ou 1: "robin" ou 2: "neumann"
    "operators_used": 1,        # 0: "basic" ou 1: "complet"
    "domain_shape": 0,          # we can choose: 0: rect - 1: img - 2: array - 3: irregular
    "final_setup" : False,       # will use f = 0 and g = y ( L - y )
    "multiplier_4_increasing_resolution": 2,
    "display" : True,
    "save" : False, 
    "image_path" : "modified_irregular05_v4.png",
    "g" : lambda x, y: 0,
}




#init_params = aff.traitement() # <---- uncomment this line to use the interface
#print(init_params)

system = SystemBuilder(**init_params)
system.build_matrixes()
system.solve_system()

# system.load_2_canva_diff_solution()       # <---- uncomment this line to display the different solutions
# system.load_2_canva_numerical_solution()  # <---- uncomment this line to display the different solutions
# system.load_2_canva_analytical_solution() # <---- uncomment this line to display the different solutions

system.load_combined_canvas()           # <---- uncomment this line to display the different solutions

plt.show()