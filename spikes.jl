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

const global rho = "rho.dat"
const global stim = "stim.dat"

const global ms = 0.001
const global sec = 1.0
const global hz = 1.0

# F = Var(SC)/Mean(SC) per window
# Spike Train, window size in ms, spike train generation duration in ms
function fano(Vector{Float64}::spikes, wid::Float64, big_t::Float64)

    # Variable Initilisation
    size_v = length(spikes)
    max_t = spikes[size_v]
    num_w = big_t/wid
    tempVector = Vector{Float64}(0,num_w)
    index = 1::Int64

    # Filter counts into tempVector
    for i in 1:num_w
        if spikes[index] < i*wid
            index++
            tempVector[i]++
        end
    end

    s_mean = mean(tempVector)
    s_var = varm(tempVector,s_mean, corrected=false)

    fanoFactor = s_var/s_mean
    return fanoFactor

end


# CV = StdDev(ISI)/Mean(ISI)
function cv(Vector{FLoat64}::spikes)

end

function interspike()

end

function queueOne()
    spikeTrain1 = [0.0]::Vector{Float64}
    spikeTrain1 = get_spike_train(35.0*hz, 1000.0*sec, 0.0*ms)
    deleteat!(spikeTrain1,1)

    #ff1 = fano(spikeTrain1)
    coeff1 = cv(spikeTrain1)

    spikeTrain2 = [0.0]::Vector{Float64}
    spikeTrain2 = get_spike_train(35.0*hz, 1000.0*sec, 5.0*ms)
    deleteat!(spikeTrain2,1)

    #ff2 = fano(spikeTrain2)
    coeff2 = cv(spikeTrain2)

end
