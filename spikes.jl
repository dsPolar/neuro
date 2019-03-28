#spikes.jl

#Q1 Calculate the Fano Factor of the SPIKE COUNT and CV of the INTER-SPIKE INTERVAL
#   for 1000 seconds of spike train with firing rate 35Hz
#   with refractory period 0 and 5ms
#   Fano Factor should use windows of width 10,50,100ms

#Q2 Calculate the Fano Factor and CV for spike train as above for rho.dat

#Q3 Calculate and plot the spike triggered average for stim.dat

#M1 Calulcate the Stimulus triggered by pairs of spikes
#   For intevals of 2,10,20,50ms calculate the average stimulus before
#   a pair of spikes seperated by that interval
#   For adjacency case and non adjacent

include("load.jl")
include("poisson.jl")
using Printf
using Statistics
using Core
using Base

const global rho = "rho.dat"
const global stim = "stim.dat"

global ms = 0.001
global sec = 1.0
global hz = 1.0

# F = Var(SC)/Mean(SC) per window
# Spike Train, window size in ms, spike train generation duration in ms
function fano(spikes::Vector{Float64}, wid::Float64, big_t::Float64)

    # Variable Initilisation
    num_w = Int(floor(big_t/wid))
    tempVector = Vector{Float64}(undef,num_w)
    index = 1::Int64

    # Filter counts into tempVector
    for i in 1:num_w
        if spikes[index] < i*wid
            index += 1
            tempVector[i] += 1
        end
    end
    #println("Calculated tempVector \n", tempVector)
    s_mean = mean(tempVector)
    s_var = varm(tempVector,s_mean)

    fanoFactor = s_var/s_mean
    return fanoFactor
end


# CV = StdDev(ISI)/Mean(ISI)
function cv(spikes::Vector{Float64})

    # Variable Initialisation
    num_s = length(spikes)
    tempVector = Vector{Float64}(undef,num_s)

    for i in 2:num_s
        tempVector[i-1] = spikes[i] - spikes[i-1]
    end
    #println("Calculated ISIs \n", tempVector)
    s_mean = mean(tempVector)
    s_std = std(tempVector, mean=s_mean)
    #println("Mean\n",s_mean,"\nStd\n",s_std)

    coeff = s_std/s_mean
    return coeff

end

function interspike()

end

function queueOne()
    spikeTrain1 = [0.0]::Vector{Float64}
    spikeTrain1 = get_spike_train(35.0*hz, 1000.0*sec, 0.0*ms)
    deleteat!(spikeTrain1,1)

    ff1 = [fano(spikeTrain1, 10.0*ms, 1000.0*ms), fano(spikeTrain1, 50.0*ms, 1000.0*ms), fano(spikeTrain1, 100.0*ms, 1000.0*ms)]
    println("Fano Factor of ST1: \n",ff1)
    coeff1 = cv(spikeTrain1)
    println("Coefficient of Variation of ST1: \n", coeff1)

    spikeTrain2 = [0.0]::Vector{Float64}
    spikeTrain2 = get_spike_train(35.0*hz, 1000.0*sec, 5.0*ms)
    deleteat!(spikeTrain2,1)

    ff2 = [fano(spikeTrain2, 10.0*ms, 1000.0*ms), fano(spikeTrain2, 50.0*ms, 1000.0*ms), fano(spikeTrain2, 100.0*ms, 1000.0*ms)]
    println("Fano Factor of ST2: \n",ff2)
    coeff2 = cv(spikeTrain2)
    println("Coefficient of Variation of ST2: \n", coeff2)

end

queueOne()
