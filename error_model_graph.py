from system_builder import SystemBuilder, plot_error_vs_h
import parameters as pr


""" 
This script is used to visualize the order of convergence of the FDM.

How to use?

The user can change the parameters of the system_builder at init_params
to check the convergence of the FDM for different problems.

The user can also change the title_name to reflect the problem being solved.
"""

title_name = "Rectangular Domain Full EDP with Robin BCs at Gamma k={pr.k}"

real_error = []
img_error = []
h = []

for i, multiplier in enumerate(range(5, 30, 5)):

    print(f"Iter {i+1} with multiplier {multiplier}")

    init_params = {
        "Lx": 10,
        "Ly": 40,
        "display": False,
        "bc_side": 1,               # 0: "dirichlet" ou 1: "robin" ou 2: "neumann"
        "operators_used": 1,        # 0: "basic" ou 1: "complet"
        "domain_shape": 0,          # we can choose: 0: rect - 1: img - 2: array
        "final_setup" : False,       # will use f = 0 and g = y ( L - y )
        "multiplier_4_increasing_resolution": multiplier,
        "display" : True,
        "save" : False,
        "image_path" : "Irregulier.png",
    }

    system = SystemBuilder(**init_params)
    system.build_matrixes()
    system.solve_system()
    
    h.append(system.metric.h)
    img_error.append(system.imaginary_error)
    real_error.append(system.real_error)   
    
plot_error_vs_h(real_error, img_error, h, title=f"{title_name}")