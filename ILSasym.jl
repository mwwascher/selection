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

function lam_2(r, p_m)
    p = (-(1-r*(1-p_m)+p_m) + sqrt((1-r*(1-p_m) + p_m)^2 - 4*(r-1)*(-p_m)))/(2*(r-1))

    p_AA = ((1-p_m)*p)/((1-p_m)*p+p_m*(1-p))
    p_Aa = 1- p_AA
    p_aA = (p_m*p)/((1-p_m)*(1-p)+p_m*p)
    p_aa = 1-p_aA

    p_00 = choose(2,0)*p_aA^0*(1-p_aA)^2
    p_01 = choose(2,1)*p_aA^1*(1-p_aA)^1
    p_02 = choose(2,2)*p_aA^2*(1-p_aA)^0

    p_10 = p_Aa*p_aa^1
    p_11 = p_AA*p_aa^1 + p_Aa*choose(1,1)*p_aA
    p_12 = p_AA*choose(1,1)*p_aA

    p_20 = p_Aa^2
    p_21 = choose(2,1)*p_AA*p_Aa
    p_22 = p_AA^2


    pmat = [p_00 p_01 p_02; p_10 p_11 p_12; p_20 p_21 p_22]
    imat = Diagonal([1,1,1])
    adj = [1; 1; 1]
    aug = hcat(pmat-imat,adj)
    p_inv = pinv(aug)
    vec = [0 0 0 1]
    pi_h = vec * p_inv

    lam_h = [0 0 0]
    lam_h = convert(Array{Float64,2}, lam_h)
    for h = 1:3
        lam_h[h] = pi_h[h]*(choose(h-1,2)*(1/p) + choose(2-(h-1),2)*(1/(1-p)))
    end
    lam = sum(lam_h)
    return lam
end

function lam_3(r, p_m)
    p = (-(1-r*(1-p_m)+p_m) + sqrt((1-r*(1-p_m) + p_m)^2 - 4*(r-1)*(-p_m)))/(2*(r-1))

    p_AA = ((1-p_m)*p)/((1-p_m)*p+p_m*(1-p))
    p_Aa = 1- p_AA
    p_aA = (p_m*p)/((1-p_m)*(1-p)+p_m*p)
    p_aa = 1-p_aA

    p_00 = choose(3,0)*p_aA^0*(1-p_aA)^3
    p_01 = choose(3,1)*p_aA^1*(1-p_aA)^2
    p_02 = choose(3,2)*p_aA^2*(1-p_aA)^1
    p_03 = choose(3,3)*p_aA^3*(1-p_aA)^0

    p_10 = p_Aa*p_aa^2
    p_11 = p_AA*p_aa^2 + p_Aa*choose(2,1)*p_aA*p_aa
    p_12 = p_AA*choose(2,1)*p_aA*p_aa^1 + p_Aa*choose(2,2)*p_aA^2
    p_13 = p_AA*choose(2,2)*p_aA^2

    p_20 = p_Aa^2*p_aa^1
    p_21 = choose(2,1)*p_AA*p_Aa*p_aa + p_Aa^2*choose(1,1)*p_aA
    p_22 = p_AA^2*p_aa^1 + choose(2,1)*p_AA*p_Aa*choose(1,1)*p_aA
    p_23 = p_AA^2*choose(1,1)*p_aA

    p_30 = p_Aa^3
    p_31 = choose(3,1)*p_AA*p_Aa^2
    p_32 = choose(3,2)*p_AA^2*p_Aa
    p_33 = p_AA^3

    pmat = [p_00 p_01 p_02 p_03; p_10 p_11 p_12 p_13; p_20 p_21 p_22 p_23; p_30 p_31 p_32 p_33]
    imat = Diagonal([1,1,1,1])
    adj = [1; 1; 1; 1]
    aug = hcat(pmat-imat,adj)
    p_inv = pinv(aug)
    vec = [0 0 0 0 1]
    pi_h = vec * p_inv

    lam_h = [0 0 0 0]
    lam_h = convert(Array{Float64,2}, lam_h)
    for h = 1:4
        lam_h[h] = pi_h[h]*(choose(h-1,2)*(1/p) + choose(3-(h-1),2)*(1/(1-p)))
    end
    lam = sum(lam_h)
    return lam
end

function lam_4(r, p_m)
    p = (-(1-r*(1-p_m)+p_m) + sqrt((1-r*(1-p_m) + p_m)^2 - 4*(r-1)*(-p_m)))/(2*(r-1))

    p_AA = ((1-p_m)*p)/((1-p_m)*p+p_m*(1-p))
    p_Aa = 1- p_AA
    p_aA = (p_m*p)/((1-p_m)*(1-p)+p_m*p)
    p_aa = 1-p_aA

    p_00 = choose(4,0)*p_aA^0*(1-p_aA)^4
    p_01 = choose(4,1)*p_aA^1*(1-p_aA)^3
    p_02 = choose(4,2)*p_aA^2*(1-p_aA)^2
    p_03 = choose(4,3)*p_aA^3*(1-p_aA)^1
    p_04 = choose(4,4)*p_aA^4*(1-p_aA)^0

    p_10 = p_Aa*p_aa^3
    p_11 = p_AA*p_aa^3 + p_Aa*choose(3,1)*p_aA*p_aa^2
    p_12 = p_AA*choose(3,1)*p_aA*p_aa^2 + p_Aa*choose(3,2)*p_aA^2*p_aa^1
    p_13 = p_AA*choose(3,2)*p_aA^2*p_aa^1 + p_Aa*choose(3,3)*p_aA^3
    p_14 = p_AA*choose(3,3)*p_aA^3

    p_20 = p_Aa^2*p_aa^2
    p_21 = choose(2,1)*p_AA*p_Aa*p_aa^2 + p_Aa^2*choose(2,1)*p_aA*p_aa
    p_22 = p_AA^2*p_aa^2 + choose(2,1)*p_AA*p_Aa*choose(2,1)*p_aA*p_aa + p_Aa^2*p_aA^2
    p_23 = p_AA^2*choose(2,1)*p_aA*p_aa + choose(2,1)*p_AA*p_Aa*p_aA^2
    p_24 = p_AA^2*p_aA^2

    p_30 = p_Aa^3*p_aa
    p_31 = choose(3,1)*p_AA*p_Aa^2*p_aa +p_Aa^3*p_aA
    p_32 = choose(3,2)*p_AA^2*p_Aa*p_aa+choose(3,1)*p_AA*p_Aa^2*p_aA
    p_33 = p_AA^3*p_aa + choose(3,2)*p_AA^2*p_Aa*p_aA
    p_34 = p_AA^3*p_aA

    p_40 = p_Aa^4
    p_41 = choose(4,1)*p_Aa^3*p_AA
    p_42 = choose(4,2)*p_Aa^2*p_AA^2
    p_43 = choose(4,3)*p_Aa^1*p_AA^3
    p_44 = p_AA^4

    pmat = [p_00 p_01 p_02 p_03 p_04; p_10 p_11 p_12 p_13 p_14; p_20 p_21 p_22 p_23 p_24; p_30 p_31 p_32 p_33 p_34; p_40 p_41 p_42 p_43 p_44]
    imat = Diagonal([1,1,1,1,1])
    adj = [1; 1; 1; 1; 1]
    aug = hcat(pmat-imat,adj)
    p_inv = pinv(aug)
    vec = [0 0 0 0 0 1]
    pi_h = vec * p_inv

    lam_h = [0 0 0 0 0]
    lam_h = convert(Array{Float64,2}, lam_h)
    for h = 1:5
        lam_h[h] = pi_h[h]*(choose(h-1,2)*(1/p) + choose(4-(h-1),2)*(1/(1-p)))
    end
    lam = sum(lam_h)
    return lam
end

function tree_sym_ISL_4(r, p_m, tau_1, tau_2, tau_3)
    coal = 0
    lam = lam_4(r, p_m)
    track_time = []
    track_coal = []
    while coal == 0
        local next_coal = dot([1 2 3 4 5 6],rand(Multinomial(1,[1/6,1/6,1/6,1/6,1/6,1/6])))
        #1 = (a,b), 2 = (a,c), 3 = (a,d), 4 = (b,c), 5 = (b,d), 6 = (c,d)
        local next_time = rand(Exponential(1/lam))
        if next_coal == 1
            min_time = tau_1
        elseif next_coal == 6
            min_time = tau_2
        else
            min_time = tau_3
        end
        if next_time > min_time
            coal = 1
            track_time = cat(track_time, [next_time], dims = (1,1))
            track_coal = cat(track_coal, [next_coal], dims = (1,1))
            #return (next_time, next_coal, min_time)
        end
    end
    tau3_1 = max(0, tau_2 - track_time[1])
    tau3_2 = max(0, tau_3 - track_time[1])
    out = tree_ISL_3(r, p_m, tau3_1, tau3_2)
    track_time = cat(track_time, [out[1]], dims = (1,1))
    track_coal = cat(track_coal, [out[2]], dims = (1,1))

    return (track_time, track_coal)
end

function tree_asym_ISL_4(r, p_m, tau_1, tau_2, tau_3)
    nextcoal = 0
    coal_time = 0
    lam2 = lam_2(r,p_m)
    lam3 = lam_3(r,p_m)
    lam4 = lam_4(r,p_m)
    time_ab = tau_1 + rand(Exponential(1/lam2))
    time_abc = tau_2 + rand(Exponential(1/lam3))
    time_abcd = tau_3 + rand(Exponential(1/lam4))
    if time_ab < tau_2
        flip = 1
        coal_time = time_ab
    elseif time_abc < tau_3
        flip = dot([1 2 3],rand(Multinomial(1,[1/3,1/3,1/3])))
        coal_time = time_abc
    else
        flip = dot([1 2 3 4 5 6],rand(Multinomial(1,[1/6,1/6,1/6,1/6,1/6,1/6])))
        coal_time = time_abcd
    end
    #1 = (a,b), 2 = (a,c), 3 = (b,c), 4 = (a,d), 5 = (b,d), 6 = (c,d)
    if flip == 1
        nextcoal = "(a,b)"
    elseif flip == 2
        nextcoal = "(a,c)"
    elseif flip == 3
        nextcoal = "(b,c)"
    elseif flip == 4
        nextcoal = "(a,d)"
    elseif flip == 5
        nextcoal = "(b,d)"
    else
        nextcoal = "(c,d)"
    end

    return (coal_time, nextcoal)
end

function tree_asym_ISL_3(r, p_m, tau_2, tau_3, time_1, coal_1)
    lam2 = lam_2(r,p_m)
    lam3 = lam_3(r,p_m)
    if coal_1 == "(a,b)"
        time_abc = max(time_1, tau_2) + rand(Exponential(1/lam2))
        time_abcd = max(time_1, tau_3) + rand(Exponential(1/lam3))
        if time_abc < tau_3
            time = time_abc
            tree = "((a,b),c),d)"
        else
            time = time_abcd
            flip = dot([1 2 3],rand(Multinomial(1,[1/3,1/3,1/3])))
            if flip == 1
                tree = "((a,b),c),d)"
            elseif flip == 2
                tree = "((a,b),d),c)"
            else
                tree = "((a,b),(c,d))"
            end
        end
    elseif coal_1 == "(a,c)"
        time_abc = time_1 + rand(Exponential(1/lam2))
        time_abcd = max(time_1, tau_3) + rand(Exponential(1/lam3))
        if time_abc < tau_3
            time = time_abc
            tree = "((a,c),b),d)"
        else
            time = time_abcd
            flip = dot([1 2 3],rand(Multinomial(1,[1/3,1/3,1/3])))
            if flip == 1
                tree = "((a,c),b),d)"
            elseif flip == 2
                tree = "((a,c),d),b)"
            else
                tree = "((a,c),(b,d))"
            end
        end
    elseif coal_1 == "(b,c)"
        time_abc = time_1 + rand(Exponential(1/lam2))
        time_abcd = max(time_1, tau_3) + rand(Exponential(1/lam3))
        if time_abc < tau_3
            time = time_abc
            tree = "((b,c),a),d)"
        else
            time = time_abcd
            flip = dot([1 2 3],rand(Multinomial(1,[1/3,1/3,1/3])))
            if flip == 1
                tree = "((b,c),a),d)"
            elseif flip == 2
                tree = "((b,c),d),a)"
            else
                tree = "((a,d),(b,c))"
            end
        end
    elseif coal_1 == "(a,d)"
        time = time_1 + rand(Exponential(1/lam3))
        flip = dot([1 2 3],rand(Multinomial(1,[1/3,1/3,1/3])))
        if flip == 1
            tree = "((a,d),b),c)"
        elseif flip == 2
            tree = "((a,d),c),b)"
        else
            tree = "((a,d),(b,c))"
        end
    elseif coal_1 == "(b,d)"
        time = time_1 + rand(Exponential(1/lam3))
        flip = dot([1 2 3],rand(Multinomial(1,[1/3,1/3,1/3])))
        if flip == 1
            tree = "((b,d),a),c)"
        elseif flip == 2
            tree = "((b,d),c),a)"
        else
            tree = "((a,c),(b,d))"
        end
    else
        time = time_1 + rand(Exponential(1/lam3))
        flip = dot([1 2 3],rand(Multinomial(1,[1/3,1/3,1/3])))
        if flip == 1
            tree = "((c,d),a),b)"
        elseif flip == 2
            tree = "((c,d),b),a)"
        else
            tree = "((a,b),(c,d))"
        end

    end
    return(time, tree)
end



function ILS_main_asym(r, p_m, tau_1, tau_2, tau_3, nreps)
    coals = []
    for i=1:nreps
        first_coal = tree_asym_ISL_4(r, p_m, tau_1, tau_2, tau_3)
        time_1 = first_coal[1]
        coal_1 = first_coal[2]
        sec_coal = tree_asym_ISL_3(r, p_m, tau_2, tau_3, time_1, coal_1)
        tree_coal = hcat(first_coal[1], first_coal[2], sec_coal[1], sec_coal[2])
        coals = cat(coals, tree_coal, dims = (1,1))
    end
    #return(coals)
    writedlm("fourtreecoal_asym.csv", coals, ',')
end
