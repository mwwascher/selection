using LinearAlgebra
using Distributions
using DelimitedFiles
using Distributed
using Random
Random.seed!(123);

cd("C:\\Users\\matth\\Documents\\JuliaWD")

function choose(n,k)
    if n < k
        return 0
    elseif n < 0
        return 0
    elseif k < 0
        return 0
    end
    res = factorial(n)/(factorial(k)*factorial(n-k))
    return res
end

##########Entries of the P matrix for the markov chain of #A lineages

# r = 2
# p_m = 0.0001
# p = (-(1-r*(1-p_m)+p_m) + sqrt((1-r*(1-p_m) + p_m)^2 - 4*(r-1)*(-p_m)))/(2*(r-1))
#
# p_AA = ((1-p_m)*p)/((1-p_m)*p+p_m*(1-p))
# p_Aa = 1- p_AA
# p_aA = (p_m*p)/((1-p_m)*(1-p)+p_m*p)
# p_aa = 1-p_aA
#
# ########4 lineages
#
# p_00 = choose(4,0)*p_aA^0*(1-p_aA)^4
# p_01 = choose(4,1)*p_aA^1*(1-p_aA)^3
# p_02 = choose(4,2)*p_aA^2*(1-p_aA)^2
# p_03 = choose(4,3)*p_aA^3*(1-p_aA)^1
# p_04 = choose(4,4)*p_aA^4*(1-p_aA)^0
#
# p_10 = p_Aa*p_aa^3
# p_11 = p_AA*p_aa^3 + p_Aa*choose(3,1)*p_aA*p_aa^2
# p_12 = p_AA*choose(3,1)*p_aA*p_aa^2 + p_Aa*choose(3,2)*p_aA^2*p_aa^1
# p_13 = p_AA*choose(3,2)*p_aA^2*p_aa^1 + p_Aa*choose(3,3)*p_aA^3
# p_14 = p_AA*choose(3,3)*p_aA^3
#
# p_20 = p_Aa^2*p_aa^2
# p_21 = choose(2,1)*p_AA*p_Aa*p_aa^2 + p_Aa^2*choose(2,1)*p_aA*p_aa
# p_22 = p_AA^2*p_aa^2 + choose(2,1)*p_AA*p_Aa*choose(2,1)*p_aA*p_aa + p_Aa^2*p_aA^2
# p_23 = p_AA^2*choose(2,1)*p_aA*p_aa + choose(2,1)*p_AA*p_Aa*p_aA^2
# p_24 = p_AA^2*p_aA^2
#
# p_30 = p_Aa^3*p_aa
# p_31 = choose(3,1)*p_AA*p_Aa^2*p_aa +p_Aa^3*p_aA
# p_32 = choose(3,2)*p_AA^2*p_Aa*p_aa+choose(3,1)*p_AA*p_Aa^2*p_aA
# p_33 = p_AA^3*p_aa + choose(3,2)*p_AA^2*p_Aa*p_aA
# p_34 = p_AA^3*p_aA
#
# p_40 = p_Aa^4
# p_41 = choose(4,1)*p_Aa^3*p_AA
# p_42 = choose(4,2)*p_Aa^2*p_AA^2
# p_43 = choose(4,3)*p_Aa^1*p_AA^3
# p_44 = p_AA^4
#
# #####Matrix stuff
# pmat = [p_00 p_01 p_02 p_03 p_04; p_10 p_11 p_12 p_13 p_14; p_20 p_21 p_22 p_23 p_24; p_30 p_31 p_32 p_33 p_34; p_40 p_41 p_42 p_43 p_44]
# imat = Diagonal([1,1,1,1,1])
# adj = [1; 1; 1; 1; 1]
# aug = hcat(pmat-imat,adj)
# p_inv = pinv(aug)
# vec = [0 0 0 0 0 1]
# pi_h = vec * p_inv
#
# lam_h = [0 0 0 0 0]
# lam_h = convert(Array{Float64,2}, lam_h)
# for h = 1:5
#     lam_h[h] = pi_h[h]*(choose(h-1,2)*(1/p) + choose(4-(h-1),2)*(1/(1-p)))
# end
# lam = sum(lam_h)
# print(lam)

#####Simulation
function prodfun_j4(N, h, p)
    prod = 1
    for j = 1:(h-1)
        prod = prod*(N*p - j + 1)
    end
    return prod
end

function prodfun_k4(N, h, p)
    prod = 1
    for k = 1:(4-h-1)
        prod = prod*(N*(1-p) - k + 1)
    end
    return prod
end

function next_h4(p_AA, p_Aa, p_aA, p_aa, h)
    if h == 0
        p_0 = choose(4,0)*p_aA^0*(1-p_aA)^4
        p_1 = choose(4,1)*p_aA^1*(1-p_aA)^3
        p_2 = choose(4,2)*p_aA^2*(1-p_aA)^2
        p_3 = choose(4,3)*p_aA^3*(1-p_aA)^1
        p_4 = choose(4,4)*p_aA^4*(1-p_aA)^0
    end
    if h == 1
        p_0 = p_Aa*p_aa^3
        p_1 = p_AA*p_aa^3 + p_Aa*choose(3,1)*p_aA*p_aa^2
        p_2 = p_AA*choose(3,1)*p_aA*p_aa^2 + p_Aa*choose(3,2)*p_aA^2*p_aa^1
        p_3 = p_AA*choose(3,2)*p_aA^2*p_aa^1 + p_Aa*choose(3,3)*p_aA^3
        p_4 = p_AA*choose(3,3)*p_aA^3
    end

    if h == 2
        p_0 = p_Aa^2*p_aa^2
        p_1 = choose(2,1)*p_AA*p_Aa*p_aa^2 + p_Aa^2*choose(2,1)*p_aA*p_aa
        p_2 = p_AA^2*p_aa^2 + choose(2,1)*p_AA*p_Aa*choose(2,1)*p_aA*p_aa + p_Aa^2*p_aA^2
        p_3 = p_AA^2*choose(2,1)*p_aA*p_aa + choose(2,1)*p_AA*p_Aa*p_aA^2
        p_4 = p_AA^2*p_aA^2
    end

    if h == 3
        p_0 = p_Aa^3*p_aa
        p_1 = choose(3,1)*p_AA*p_Aa^2*p_aa +p_Aa^3*p_aA
        p_2 = choose(3,2)*p_AA^2*p_Aa*p_aa+choose(3,1)*p_AA*p_Aa^2*p_aA
        p_3 = p_AA^3*p_aa + choose(3,2)*p_AA^2*p_Aa*p_aA
        p_4 = p_AA^3*p_aA
    end

    if h == 4
        p_0 = p_Aa^4
        p_1 = choose(4,1)*p_Aa^3*p_AA
        p_2 = choose(4,2)*p_Aa^2*p_AA^2
        p_3 = choose(4,3)*p_Aa^1*p_AA^3
        p_4 = p_AA^4
    end
    index = rand(Multinomial(1,[p_0, p_1, p_2, p_3, p_4]))
    h_next = dot([0 1 2 3 4],index)
    return h_next
end

function pfun_i4(T, N, r, p_m, p_init)
    p_i = zeros(T)
    p_i[1] = p_init

    for i = 1:(T-1)
        X = rand(Binomial(N, (r*p_i[i])/(r*p_i[i] + (1 - p_i[i]))))
        Z_1 = rand(Binomial(X, 1-p_m))
        Z_2 = rand(Binomial(N-X, p_m))
        Y = Z_1 + Z_2
        p_i[i+1] = Y/N
    end
    return p_i
end

function timefun_c4(T, N, r, p_m, p_init)
    p_i = pfun_i4(T, N, r, p_m, p_init)
    flip = 0
    time_c = T
    h = rand(Binomial(4,p_i[T]))
    track_h_local = []

    while flip == 0 && time_c > 0
        prod_1 = prodfun_j4(N, h, p_i[time_c])
        prod_2 = prodfun_k4(N, h, p_i[time_c])

        p_flip = (choose(h,2)*prod_1)/((N*p_i[time_c])^h) + (choose(4-h,2)*prod_2)/((N*(1-p_i[time_c]))^(4-h))
        flip = rand(Binomial(1,p_flip))

        p_AA = ((1-p_m)*p_i[time_c])/((1-p_m)*p_i[time_c]+p_m*(1-p_i[time_c]))
        p_Aa = 1- p_AA
        p_aA = (p_m*p_i[time_c])/((1-p_m)*(1-p_i[time_c])+p_m*p_i[time_c])
        p_aa = 1-p_aA

        h = next_h4(p_AA, p_Aa, p_aA, p_aa, h)
        time_c = time_c - 1
        track_h_local = cat(track_h_local,[h], dims = (1,1))
    end
    return(T-time_c, track_h_local)
end

function time_sim4(r, p_m, N, nreps)
    T = 2*N
    p_init = .5

    time_hist = []
    track_h = []

    for l = 1:nreps
        out = timefun_c4(T, N, r, p_m, p_init)
        time_hist = cat(time_hist,[out[1]],dims = (1,1))
        track_h = cat(track_h,out[2], dims = (1,1))
        #print(time_c)
        print(l)
    end
    writedlm("fourtreetimes.csv", time_hist, ',')
    writedlm("fourtreeh.csv", track_h, ',')
    return nothing
end

#syntax is time_sim4(r, p_m, N, nreps)
time_sim4(1.0002, .0001, 500000, 10)
