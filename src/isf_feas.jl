#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Version 0.1, 2025.
#
# If you found this piece of software useful, please cite the paper:
# J.-P. Dussault, J.Ch. Gilbert, B. Plaquevent-Jourdain,
# "Primal and Dual Approaches for the Chamber Enumeration
# of real hyperplane arrangements", 2025.
#
# Authors:
# - Jean-Pierre Dussault (Univ. of Sherbrooke, Canada),
# - Jean Charles Gilbert (INRIA, France),
# - Baptiste Plaquevent-Jourdain (INRIA & Univ. of Sherbrooke, Canada).
#
# Copyright 2025, INRIA (France) and Univ. of Sherbrooke (Canada).
#
# ISF is distributed under the terms of the Q Public License version
# 1.0.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Q Public
# License version 1.0 for more details.
#
# You should have received a copy of the Q Public License version 1.0
# along with this program. If not, see
# <https://doc.qt.io/archives/3.3/license.html>.
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

### 
# Depending on the case, the computations are a bit different. 
# The general framework is still the same. 

"""
`(d|x), lambda = isf_feas!(VS, sv, d0, info, options, values)`

Returns the direction, the multiplier (if used) and modifies info. 
Verification of the inequalities system through linear optimization.
A version covers the mixed (HnH) algorithm, this version covers the rest.

First, the elements of the linear system are created, then the 
properties of Gurobi are set, and finally the result is checked. 
Gurobi has (apparently) the LP convention that the dual 
variables have the same sign as contraints: for this 
reason, the constraints have been put in the other sense.

The HnH version is done separately for convenience.
"""
function isf_feas!(VS::Matrix, sv::Vector, rhs, info::Info, options::Options, values::Values)

    # dimensions
    if options.symmetry
        n, km = size(VS)
        b = ones(km+2)                                      # homogeneity so constraints are s_i v_i' * d >= 1
        b[km+1] = 0                                         # no right-and side for this constraint
        b[km+2] = -1                                        # (arbitrary) bound of -1
        # LOP matrix
        A = [[VS' ; sv' ; zeros(1,n)] [zeros(km,1);1;1]]
    else
        km = size(VS,2)
        n = size(VS,1)-1
        b       = zeros(km+2)
        b[1:km] = VS[n+1,1:km]                              # there is no homogeneity so one cannot put ones
        b[km+1] = +rhs                                      # takes into account the appropriate inequality, so later the optimal value is compared with 0
        b[km+2] = -1                                        # (arbitrary) bound of -1
        A = [[VS[1:n,:]' ; sv' ; zeros(1,n)] ones(km+2,1)]  # the linear problem cannot be the same due to the lack of homogeneity, 't' must intervene
    end

    x = Inf * ones(n)                                       # by default infinity, if there is a direction found it will update and be returned
    λ = zeros(km)                                           # by default zero, if used later it will update and be returned
    k = km+1
    
    # LOP data
    c       = zeros(n+1)
    c[n+1]  = 1

    info.nb_losolve += 1

    model = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
    # model = Model(GLPK.Optimizer)                         # for GLPK
    set_silent(model)                                       # avoid printings
    set_optimizer_attribute(model, "OutPutFlag", 0)         # ?
    set_optimizer_attribute(model, "Method", -1)            # dual simplex
    # set_optimizer_attribute(model, "Dual", 1)             # for GLPK
    @variable(model, xt[1:n+1])                             # dimension + 1 for the bounding variable
    @constraint(model, cons, A * xt >= b)                   # linear constraint
    @objective(model, Min, c' * xt)                         # cost = last variable, the bounding variable

    optimize!(model)
    
    if objective_value(model) < 0                           # feasible system have an optimal value (-1) < 0
        info.nb_feaslop += 1
        x = JuMP.value.(xt)[1:n]
    else                                                    # if infeasible one takes the dual variables (and x stays at Inf * ones(n))
        info.nb_infeaslop += 1
        λ = dual(cons)[1:k]
    end

    info.flag = values.success

    return x, λ
end

#-----------------------------------------------------------------------

"""
`d, λ =   isf_feas_HnH!(VS, sv, d0, info, options, values)`

Returns the direction, the multiplier (if used) and modifies info. 
Verification of the inequalities system through linear optimization.
Some slight modifications arise to treat the HnH case.

First, the elements of the linear system are created, 
then the properties of Gurobi are set, and finally
the result is checked. 
Gurobi has (apparently) the LP convention that the dual 
variables have the same sign as contraints: for this 
reason, the constraints have been put in the other sense.
"""
function isf_feas_HnH!(VS::Matrix, sv::Vector, info::Info, options::Options, values::Values)

    km = size(VS, 2)
    n  = size(VS, 1)                                        # already the right dimension by definition of the isf_feas_HnH!(...) calls

    b = ones(km+2)                                          # homogeneity so constraints are s_i v_i' * d >= 1
    b[km+1] = 0                                             # no right-and side for this constraint
    b[km+2] = -1                                            # (arbitrary) bound of -1
    A = [[VS' ; sv' ; zeros(1,n)] [zeros(km,1);1;1]]

    x = Inf * ones(n)                                       # by default infinity, if there is a direction found it will update and be returned
    λ = zeros(km)                                           # by default zero, if used later it will update and be returned
    k = km+1
    
    # LOP data
    c       = zeros(n+1)
    c[n+1]  = 1

    info.nb_losolve += 1

    model = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
    # model = Model(GLPK.Optimizer)                         # for GLPK
    set_silent(model)                                       # avoid printings
    set_optimizer_attribute(model, "OutPutFlag", 0)         # ?
    set_optimizer_attribute(model, "Method", -1)            # dual simplex
    # set_optimizer_attribute(model, "Dual", 1)             # for GLPK
    @variable(model, xt[1:n+1])                             # dimension + 1 for the bounding variable
    @constraint(model, cons, A * xt >= b)                   # linear constraint
    @objective(model, Min, c' * xt)                         # cost = last variable, the bounding variable

    optimize!(model)
    
    if objective_value(model) < 0                           # feasible system, value should be -1 < 0
        info.nb_feaslop += 1
        x = JuMP.value.(xt)[1:n]
    else                                                    # if infeasible one takes the dual variables (and x stays at Inf * ones(n))
        info.nb_infeaslop += 1
        λ = dual(cons)[1:k]
    end

    info.flag = values.success

    return x, λ
end