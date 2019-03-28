#spikes.jl

#Q1 Calculate the Fano Factor of the SPIKE COUNT and CV of the INTER-SPIKE INTERVAL
#   for 1000 dSeconds of spike train with firing rate 35dHz
#   with refractory period 0 and 5dMs
#   Fano Factor should use windows of width 10,50,100dMs

#Q2 Calculate the Fano Factor and CV for spike train as above for rho.dat
#   Requires mapping rho.dat to spike times, 2dMs sampling, 500dHz

#Q3 Calculate and plot the spike triggered average for stim.dat

#M1 Calulcate the Stimulus triggered by pairs of spikes
#   For intevals of 2,10,20,50dMs calculate the average stimulus before
#   a pair of spikes seperated by that interval
#   For adjacency case and non adjacent

include("load.jl")
include("poisson.jl")
using Printf
using Statistics
using PyPlot

const global rhodat = "rho.dat"
const global stimdat = "stim.dat"

const global dMs = 0.001
const global dSec = 1.0
const global dMin = 60.0
const global dHz = 1.0

# F = Var(SC)/Mean(SC) per window
# Spike Train, window size in dMs, spike train generation duration in dMs
function fano(spikes::Vector{Float64}, wid::Float64, big_t::Float64)

    # Variable Initilisation
    num_w = Int(floor(big_t/wid))
    #println("Big_t",big_t)
    #println("wid", wid)
    #println("Num_w",num_w)
    tempVector = zeros(num_w)
    win_num = 1::Int64

    for i in 1:length(spikes)
        if spikes[i] < (win_num * wid)
            tempVector[win_num] += 1
        else
            win_num += 1
            tempVector[win_num] += 1
        end
    end


    for z in 1:length(tempVector)
        if tempVector[z] > 1000.0
            #println("Insanity at index: \n",z)
            #println("Insanity value: \n", tempVector[z])
            #if z+1 > length(tempVector)
            #    println("At end of vector")
            #end
            tempVector[z] = 0.0
        end
    end

    s_mean = mean(tempVector)
    s_var = varm(tempVector,s_mean)

    fanoFactor = s_var/s_mean
    return fanoFactor
end


# CV = StdDev(ISI)/Mean(ISI)
function cv(spikes::Vector{Float64})

    # Variable Initialisation
    num_s = length(spikes)
    tempVector = zeros(num_s)

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


function convertRho()
    spikes = loadRho()
    count = 0::Int64

    for c in 1:length(spikes)
        if spikes[c] == 1
            count += 1
        end
    end

    spikeTrain = zeros(count)
    index = 1::Int64

    for i in 1:length(spikes)
        if spikes[i] == 1
            spikeTrain[index] = spikes[i]*i*2*dMs
            index += 1
        end
    end
    return spikeTrain
end

function loadRho()
    rho = load_data(rhodat,Int64)
    return rho
end

function loadStim()
    stimulus = load_data(stimdat,Float64)
    return stimulus
end

function spikeTriggered()
    rho = [0]::Vector{Int64}
    rho = loadRho()
    stimulus = [0.0]::Vector{Float64}
    stimulus = loadStim()



    count = 0::Int64
    for i in 1:length(rho)
        if rho[i] == 1
            count += 1
        end
    end

    averages = zeros(50)

    for i in 1:length(rho)
        if rho[i] == 1
            index_s = i
            for bt in 1:(50)
                if i-bt > 0
                    averages[51-bt] += stimulus[i-bt]
                end
            end
        end
    end

    newAverages = zeros(100)

    for i in 1:length(averages)
        averages[i] = averages[i] / count
    end

    return averages
end


function queueOne()
    spikeTrain1 = [0.0]::Vector{Float64}
    spikeTrain1 = get_spike_train(35.0*dHz, 1000.0*dSec, 0.0*dMs)
    deleteat!(spikeTrain1,1)

    ff1 = [fano(spikeTrain1, 10.0*dMs, 1000.0*dSec), fano(spikeTrain1, 50.0*dMs, 1000.0*dSec), fano(spikeTrain1, 100.0*dMs, 1000.0*dSec)]
    println("Fano Factor of ST1: \n",ff1)
    coeff1 = cv(spikeTrain1)
    println("Coefficient of Variation of ST1: \n", coeff1)

    spikeTrain2 = [0.0]::Vector{Float64}
    spikeTrain2 = get_spike_train(35.0*dHz, 1000.0*dSec, 5.0*dMs)
    deleteat!(spikeTrain2,1)

    ff2 = [fano(spikeTrain2, 10.0*dMs, 1000.0*dSec), fano(spikeTrain2, 50.0*dMs, 1000.0*dSec), fano(spikeTrain2, 100.0*dMs, 1000.0*dSec)]
    println("Fano Factor of ST2: \n",ff2)
    coeff2 = cv(spikeTrain2)
    println("Coefficient of Variation of ST2: \n", coeff2)
end

function queueTwo()
    spikeTrain = [0.0]::Vector{Float64}
    spikeTrain = convertRho()
    deleteat!(spikeTrain,1)

    ff = [fano(spikeTrain, 10.0*dMs, 20*dMin*dSec), fano(spikeTrain, 50.0*dMs, 20*60*dSec), fano(spikeTrain, 100.0*dMs, 20*60*dSec)]
    println("Fano Factor of rho.dat: \n",ff)

    coeff = cv(spikeTrain)
    println("Coefficient of Variation of rho.dat: \n",coeff)
end

function queueThree()

    averages = spikeTriggered()
    #println("STA Averages: \n", averages)
    stimulus = [0.0]::Vector{Float64}
    stimulus = loadStim()

    plot(averages)
    plt.title("Spike Triggered Average plot for Rho.dat and Stim.dat\n")
    plt.xlabel("Sample window size/samples (1 per 2ms)")
    plt.ylabel("Stimulus Value")
    savefig("sta.svg")
    println("File saved to sta.svg")

end


queueOne()
queueTwo()
queueThree()
