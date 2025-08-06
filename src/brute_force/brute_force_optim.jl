### 
# d, val =   brute_force_optim_H(V)
# d, val =  brute_force_optim_nH(V)

# Returns a direction d sch that V' * d > 0, and [] otherwise,
# done by checking the value of the linear optimization problem
#
# { min t 
# { V' *d + t*e >= 0
# { t >= -1
# 
# is negative (e = [1,...,1]) 

# Multiple (very similar) variants are used for the h, nH, HnH cases.
# The details change a little bit

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Version 0.1, May, 2023.
#
# If you found this piece of software useful, please cite the paper:
# J.-P. Dussault, J.Ch. Gilbert, B. Plaquevent-Jourdain,
# '?????', 2023.
#
# Authors:
# - Jean-Pierre Dussault (Univ. of Sherbrooke, Canada),
# - Jean Charles Gilbert (INRIA, France),
# - Baptiste Plaquevent-Jourdain (INRIA & Univ. of Sherbrooke, Canada).
#
# Copyright 2023, INRIA (France) and Univ. of Sherbrooke (Canada).
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
# <http://doc.trolltech.com/3.0/license.html>.
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

function brute_force_optim_H(V)

    n, p = size(V)

    # LOP data
    c      = zeros(n+1)
    c[n+1] = 1
    A      = -[[V' ; zeros(1,n)] ones(p+1,1)]
    rhs      = zeros(p+1)
    rhs[p+1] = 1

    model = Model(Gurobi.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "Method", 1);
    @variable(model, dt[1:n+1])
    @constraint(model, cons, A * dt <= rhs)
    @objective(model, Min, c' * dt)

    optimize!(model)

    d = JuMP.value.(dt)[1:n]
    t = JuMP.value.(dt)[n+1]

    # flag shenanigans

    return d, t
end

function brute_force_optim_nH(Vt)
    n = size(Vt, 1) - 1
    p = size(Vt, 2)
    
    V = Vt[1:n,:]
    T = Vt[n+1,:]

    c      = zeros(n+1)
    c[n+1] = 1
    A      = -[[V' ; zeros(1,n)] ones(p+1,1)]
    rhs      = zeros(p+1)
    rhs[1:p] = -T
    rhs[p+1] = 1

    model = Model(Gurobi.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "Method", 1);
    @variable(model, dt[1:n+1])
    @constraint(model, cons, A * dt <= rhs)
    @objective(model, Min, c' * dt)

    optimize!(model)

    d = JuMP.value.(dt)[1:n]
    t = JuMP.value.(dt)[n+1]

    # flag shenanigans

    return d, t
end
