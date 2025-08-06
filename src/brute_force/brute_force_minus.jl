### 
# b = brute_force_minus(b)

# Updates the binary vector:
# decreases '1' from b binarywise; note that the first bit is the
# less significant one. It returns [] if b = 0

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

function brute_force_minus(binary)

    # returns 0 if b == 0
    if all(binary .== 0)
        binary = []
        return binary
    end

    # decrease 1

    n = length(binary)
    for i in 1:n
        if binary[i] == 0
            binary[i] = 1
        else
            binary[i] = 0
            break
        end
    end

    return binary
end