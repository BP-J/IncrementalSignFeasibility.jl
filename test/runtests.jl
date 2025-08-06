using IncrementalSignFeasibility
using Test

# short tests for the package
# more detailed tests exist (to verify the paper's results)

@testset "IncrementalSignFeasibility.jl" begin
    @testset "Correctness tests" begin
        rand_12_6 = readdlm("affine_data/data_aff_R_12_6.txt")
        twod_6 = readdlm("affine_data/data_aff_degen2d_6.txt")
        srand_2 = readdlm("affine_data/data_aff_srand_8_20_2.txt")
        perm_6 = readdlm("affine_data/data_aff_perm_6.txt")
        ratio_5_7 = readdlm("affine_data/data_aff_ratio_5_7.txt")

        options_7 = Options(7, 0, false, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, false)  # versions with linear optimization and Gurobi seemed to cause an error in the test (not the algorithm)
        info_rand_12_6 = isf(rand_12_6, options_7)
        info_twod_6 = isf(twod_6, options_7)
        info_srand_2 = isf(srand_2, options_7)
        info_perm_6 = isf(perm_6, options_7)
        info_ratio_5_7 = isf(ratio_5_7, options_7)
        @test [info_rand_12_6.ns, info_twod_6.ns, info_srand_2.ns, info_perm_6.ns, info_ratio_5_7.ns]==[2510, 2176, 36225, 5040, 15136]
    end
end
