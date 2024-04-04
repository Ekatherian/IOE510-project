import gurobipy as gp
from gurobipy import GRB

# Initialize the model
model = gp.Model("Power System and Wildfire Risk Management")

#################################################
####                   Sets                  ####
#################################################
# set of subareas
A_e = {'Area1', 'Area2'}

# set of days in a week
W = range(7)
# set of hours in a day
H = range(24)

# Loads
D_wh = {'Load1':100, 'Load2':150}
# Generators
G = {'Gen1', 'Gen2'} 
# Transmission lines
L = {'Line1', 'Line2'} 
# Buses 
B = {'Bus1', 'Bus2'}  
# electrical components, union of G,L,B and D_wh
E = G.union(L, B, D_wh)

#################################################
####              Parameters                 ####
#################################################

alpha = 0.5  # Trade-off parameter

# Base risk of wildfire ignition in each area
R_j = {'Area1': 1.2, 'Area2': 1.5}

# Relative risk factor for electrical components in areas at each hour
kappa_ejh = {('Gen1', 'Area1', 1): 1.1, ('Line1', 'Area2', 1): 1.2}

# Meteorological factor at each hour
gamma_h = {1: 1.2, 2: 1.1}

# Risk calculations
R_ejh = {(e, j, h): kappa_ejh.get((e, j, h), 1.0) * gamma_h.get(h, 1.0) * R_j[j] # 
         for e in E for j in A_e for h in H}
# The .get() method is used to retrieve the value from the dictionary kappa_ejh for the key (e, j, h). 
# If the key doesn’t exist, it returns 1.0 as a default value. 

R_eh = {(e,h) : sum(R_ejh[e, j, h] for j in A_e if (e, j, h) in R_ejh)
        for e in E for h in H}

# Total wildfire risk calculation
R_Fire_h = {h: sum(x_dh[e, h] * R_ejh.get((e, j, h), 0) for e in D_wh for j in A_e)
            + sum(z_gh[e,h] * R_ejh.get((e,j,h),0) for e in G for j in A_e)
            + sum(z_lh[e,h] * R_ejh((e,j,h),0) for e in L for j in A_e)
            + sum(z_ih[e,h] * R_ejh((e,j,h),0) for e in B for j in A_e)
            for h in H}

# w_dh represents the weighting factor for each load and hour, initializing to 1
w_dh = {(d, h): 1.0 for d in D_wh for h in H} 

# Typical load served during standard operational phases for each load and hour
D_dh = {(d, h): D_wh[d] for d in D_wh for h in H}

# Calculate the total load successfully distributed across the network within hour h
D_Tot_h = {h: sum(x_dh[d, h] * w_dh[d, h] * D_dh[d, h] for d in D_wh) for h in H}

# Define generator min/max power output
P_min, P_max = 10, 100  # Example values

T_l = {l: 50 for l in L}  # Example thermal limits
b_l = {l: 0.01 for l in L}  # Example susceptance values
theta_max, theta_min = 30, -30  # Example angle limits



#################################################
####           Decision Variables            ####
#################################################

# indicate whether load l will remain energized during hour h
z_lh = model.addVars(L, H, vtype=GRB.BINARY, name="z_lh")

# indicate whether bus i will remain energized during hour h
z_ih = model.addVars(B, H, vtype=GRB.BINARY, name="z_ih")

# indicate whether generator l will remain energized during hour h
z_gh = model.addVars(G, H, vtype=GRB.BINARY, name="z_gh")

# represent the fraction of the load that that is de-energized during hour h
x_dh = model.addVars(A_e, H, vtype=GRB.CONTINUOUS, name="x_dh", lb=0, ub=1)

# represent the power generation of generators in hour h
P_gh = model.addVars(G, H, lb=0, name="P_gh")

# represent the power flows on transmission lines l in L from bus j to bus i in hour h
P_l_ij_h = model.addVars(L, B, B, H, lb=-GRB.INFINITY, name="P_l_ij_h")

# represent the voltage angles.
theta_ih = model.addVars(B, H, lb=-GRB.INFINITY, ub=GRB.INFINITY, name="theta_ih")


#################################################
####           Objective Function            ####
#################################################

model.setObjective(alpha * sum(D_Tot_h[h] for h in H) - (1 - alpha) * sum(R_Fire_h[h] for h in H), GRB.MAXIMIZE)



#################################################
####               Constraints               ####
#################################################


# Energization constraints
for i in B:
    for h in H:
        # Loads must be less than or equal to bus energization status
        model.addConstrs((x_dh[d, h] <= z_ih[i, h] for d in D_wh if i in B), name=f"load_energization_{i}_{h}")
        # Generator status must be less than or equal to bus energization status
        model.addConstrs((z_gh[g, h] <= z_ih[i, h] for g in G if i in B), name=f"gen_energization_{i}_{h}")
        # Line energization
        model.addConstrs((z_lh[l, h] <= z_ih[i, h] for l in L if i in B), name=f"line_energization_{i}_{h}")

# Generation constraints
for g in G:
    for h in H:
        model.addConstr(P_min * z_gh[g, h] <= P_gh[g, h], name=f"Pmin_{g}_{h}")
        model.addConstr(P_gh[g, h] <= P_max * z_gh[g, h], name=f"Pmax_{g}_{h}")

# Power flow and node balance constraints
for i in B:
    for h in H:
        model.addConstr(
            sum(P_gh[g, h] for g in G if i in B) +
            sum(P_l_ij_h[l, i, j, h] for l in L for j in B if i in B) -
            sum(x_dh[d, h] * D_wh[d] for d in D_wh if i in B) == 0,
            name=f"node_balance_{i}_{h}"
        )

        for l in L:
            for j in B:
                if i != j:  # Only add flow constraints for lines between different buses
                    model.addConstr(
                        -T_l[l] * z_lh[l, h] <= P_l_ij_h[l, i, j, h],
                        name=f"flow_limit_lower_{l}_{i}_{j}_{h}"
                    )
                    model.addConstr(
                        P_l_ij_h[l, i, j, h] <= T_l[l] * z_lh[l, h],
                        name=f"flow_limit_upper_{l}_{i}_{j}_{h}"
                    )
                    model.addConstr(
                        P_l_ij_h[l, i, j, h] <= -b_l[l] * (theta_ih[i, h] - theta_ih[j, h] + theta_max * (1 - z_lh[l, h])),
                        name=f"angle_limit_upper_{l}_{i}_{j}_{h}"
                    )
                    model.addConstr(
                        P_l_ij_h[l, i, j, h] >= -b_l[l] * (theta_ih[i, h] - theta_ih[j, h] + theta_min * (1 - z_lh[l, h])),
                        name=f"angle_limit_lower_{l}_{i}_{j}_{h}"
                    )

# Solve the model
model.optimize()


# Check if the optimal solution is found
if model.status == GRB.OPTIMAL:
    print("Optimal Solution:")

    # Print operational status for generators, lines, and buses
    for h in H:
        print(f"\nHour {h}:")
        for g in G:
            print(f"  Generator {g}: Energized = {'Yes' if z_gh[g, h].X > 0.5 else 'No'}, Power Generated = {P_gh[g, h].X:.2f} MW")
        for l in L:
            print(f"  Line {l}: Energized = {'Yes' if z_lh[l, h].X > 0.5 else 'No'}")
        for b in B:
            print(f"  Bus {b}: Energized = {'Yes' if z_ih[b, h].X > 0.5 else 'No'}")

        # Print load shedding details
        for d in D_wh:
            print(f"  Load {d} Shedding: {x_dh[d, h].X * 100:.2f}% of {D_wh[d]} MW")

        # Print power flows on transmission lines if needed
        for l in L:
            for i in B:
                for j in B:
                    if i != j and P_l_ij_h[l, i, j, h].X != 0:
                        print(f"  Power Flow on Line {l} from Bus {i} to Bus {j}: {P_l_ij_h[l, i, j, h].X:.2f} MW")

        # Print voltage angles at buses if needed
        for b in B:
            print(f"  Voltage Angle at Bus {b}: {theta_ih[b, h].X:.2f} degrees")
else:
    print("No optimal solution found.")
