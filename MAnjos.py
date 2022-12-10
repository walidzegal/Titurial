"""
Created on Wed Nov 30 16:36:28 2022
@author: ZW
"""
#TRANSPORTATION MODEL formulation
from pyomo.environ import *
#dexlare model
model= AbstractModel()

# sets
model.SUPPLIES= Set()#set of supply nodes
model.T= Set()#SET OF DEMAND NODES

#PARAMS
model.s = Param(model.SUPPLIES) # the supply available at each node
model.d = Param(model.T) # the demand required at t
model.c = Param(model.SUPPLIES, model.T) # the cost of sipping from supply node eo demand node
model.R = Param(model.T) # the reserve requiment at t
model.cu = Param(model.SUPPLIES)#startup cost
model.Ru=Param(model.SUPPLIES) #maximum ramp-up rate 
model.Su=Param(model.SUPPLIES) # maximum startup rate
model.Rd=Param(model.SUPPLIES) # maximum ramp-down rate
model.Sd = Param(model.SUPPLIES)# maximum shutdown rate
model.Tu=Param(model.SUPPLIES)#min number of periods on
model.Td=Param(model.SUPPLIES)#min number of periods off
model.pUp= Param(model.SUPPLIES)#power upper bound
model.plow = Param(model.SUPPLIES)# power lower bound

#decision variable x[i,j] = amount shipped from plant i at time t
model.p = Var(model.SUPPLIES,model.T,within=NonNegativeReals)
model.pup  = Var(model.SUPPLIES,model.T,within=NonNegativeReals)
model.y = Var(model.SUPPLIES,model.T,within=Binary)# y[j,t]=1 if j start up at t, 0 otherwise
model.v = Var(model.SUPPLIES,model.T,within=Binary)# v[j,t]=1 if j is on at t, 0 otherwise
model.z = Var(model.SUPPLIES,model.T,within=Binary)# z[j,t]=1 if j is shuts down at beginning of t, 0 otherwise
# before horizon
# =============================================================================
# for j in model.SUPPLIES :
#     model.p[j,0].setlb(0)
#     model.p[j,0].setun(0)
#     model.p[1,0].setlb(120)
#     model.p[1,0].setub(120)
# =============================================================================
#objective function : minimize totalcost to meet demand
def obejective_rule(model):
    return sum(model.c[i,t]*model.p[i,t] +model.cu[i]*model.y[i,t]for i in model.SUPPLIES for t in model.T-{0})
model.minCost= Objective(rule=obejective_rule,sense=minimize)

#demand constraints :
def demand_rule(model, t):
    return (sum(model.p[i,t]for i in model.SUPPLIES)==model.d[t])
model.demandConstraints =Constraint(model.T,rule=demand_rule)

 # #reserve constraints :
def reserve_rule(model,t) :
     return(sum(model.pup[j,t] for j in model.SUPPLIES)>=model.d[t]+model.R[t])
model.reserveConstraints = Constraint(model.T-{1},rule=reserve_rule)

# # logical coherence constraints :
def logical_rule(model,j,t) :
     return(model.v[j,t-1]-model.v[j,t]+model.y[j,t]-model.z[j,t]==0 )
model.logical_constraints= Constraint(model.SUPPLIES,model.T-{1},rule=logical_rule)


## Generation limits ramping, startup and downtime rate
def ramping_startup_rule(model,j,t):
    return(model.p[j,t]-model.p[j,t-1]<=model.Ru[j]*model.v[j,t-1]+model.Su[j]*model.y[j,t] )
model.ramping_startup_constraints = Constraint(model.SUPPLIES,model.T-{1},rule=ramping_startup_rule)


# ramping a constraints with shutdown rate
    
def ramping_shutdown_rule(model,j,t):
    
    return(model.p[j,t-1]-model.p[j,t]<= model.Rd[j]*model.v[j,t]+model.Sd[j]*model.z[j,t])
model.ramping_shutdown_constraints = Constraint(model.SUPPLIES,model.T-{1},rule=ramping_shutdown_rule)

    
#uptime constraint
def uptime_rule(model,j,t):
    k1 = max(t-model.Tu[j]+1,1)
    return(sum(model.y[j,k] for k in range(k1,t))<=model.v[j,t])

model.uptime_constraints = Constraint(model.SUPPLIES,model.T-{0},rule=uptime_rule)
    
#shut down time constraint
def shut_down_time_rule(model,j,t):
    k1 = max(t-model.Td[j]+1,1)
    return(model.v[j,t]+sum(model.z[j,k] for k in range(k1,t))<=1)

model.shuttime_constraints = Constraint(model.SUPPLIES,model.T,rule=shut_down_time_rule)


## generation limits
def genaeration_limits1_rule(model,j,t):
    return(model.plow[j]*model.v[j,t]<=model.p[j,t])
model.generation_limits_constraints= Constraint(model.SUPPLIES,model.T,rule=genaeration_limits1_rule)
## generation limits
def genaeration_limits2_rule(model,j,t):
    return(model.p[j,t]<=model.pUp[j]*model.v[j,t])
model.generation_limits2_constraints= Constraint(model.SUPPLIES,model.T,rule=genaeration_limits2_rule)

#generation startuo limts
def genaeration_startuplimits_rule(model,j,t):
    return(model.pup[j,t]<=model.p[j,t-1]+model.Ru[j]*model.v[j,t-1]+model.Su[j]*model.y[j,t])
model.generation_startuplimits_constraints= Constraint(model.SUPPLIES,model.T-{1},rule=genaeration_startuplimits_rule)

#generation shutdown limts
def genaeration_shutdownlimits_rule(model,j,t):
    return(model.pup[j,t]<=model.pUp[j]*(model.v[j,t]-model.z[j,t+1])+model.Sd[j]*model.z[j,t+1])
model.generation_shutdownlimits_constraints= Constraint(model.SUPPLIES,model.T-{6},rule=genaeration_shutdownlimits_rule)

# # creating a model instance (combining the abstract model with a specific data file
data= DataPortal()
data.load(filename="MAnjos.dat", model=model)
instance = model.create_instance(data)

# # Solving the instance
optimiser = SolverFactory('glpk')
optimiser.solve(instance)
instance.display()
print('optimal value :', value(instance.minCost))
