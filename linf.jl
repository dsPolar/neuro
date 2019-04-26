#CW3

#Q1

# Simulate an integrate and fire model with parameters below for one second
# tau_m = 10ms       Membrane time constant
# E_L = V_r = -70mV  Leak Potential, Reset Voltage
# V_t = -40mV        Threshold voltage
# R_m = 10M          Membrane Resistance (receprocal of conductance)
# I_e = 3.1nA

# Use Euler's method with time step delta_t = 1ms
# Plot the voltage as a function of time
# No Refractory period
# Once voltage exceeds V_t set voltage to V_r


using Printf
using PyPlot
using Distributions

const global tau_m = 0.01
const global E_L = -0.07
const global V_r = -0.07
const global V_t = -0.04
const global R_m = 10000000
const global C_m = 1/R_m
const global I_e = 0.0000000031
const global delta_t = 0.001



#Update a voltage amount for time step
function voltageUpdate(voltage::Float64)
    #Taken from conor's notes
    # tau_m dV/dt = E_L - V + R_m*I_e
    update = ((E_L - voltage) + (R_m*I_e))
    update = update * tau_m

    newVoltage = voltage + update

    return newVoltage
end

#Check if voltage is above V_t
function checkReset(voltage::Float64)
    if voltage > V_t
        return true
    else
        return false
    end
end

#Use pyplot to plot voltages vector and save to file
function plotVoltage(voltages::Vector{Float64}, path::String)
    plot(voltages)
    plt.title("Modelled Leaky Integrate and Fire Neuron\n")
    plt.xlabel("time/ms")
    plt.ylabel("Voltage/V")
    savefig(path)
    println("File saved to ",path)
end



#Model a single linf neuron and plot voltage against time
function q1()
    voltage = zeros(1001)
    voltage[1] = V_r
    # 2:1001 to account for t[0] stored at voltage[1] due to Julia
    for t in 2:1001
        voltage[t] = voltageUpdate(voltage[t-1])
        if checkReset(voltage[t])
            voltage[t] = V_r
        end

    end
    plotVoltage(voltage, "linf.jpg")
end









#Question 2
# Simulate two neurons with synaptic connection between each other
# First projects to second
# Second projects to first

const global tau_m2 = 0.02
const global V_r2 = -0.08
const global V_t2 = -0.054
const global R_mI_e = 0.018
const global R_mg_s = 0.15
const global P = 0.5
const global tau_s = 0.010

#Simulate two cases
# a) assuming synapses are excitatory with E_s = 0
# b) assuming synapses are inhibitory with E_s = -0.08

#Set the inital membrane potentials of neurons to different random values V_r:V_t
#Simulate 1s of activity
#Plot each case plot neurons together differentiated by colour


function checkReset2(voltage::Float64)
    if voltage > V_t2
        return true
    else
        return false
    end
end

function updateVoltages2(v1::Float64, v2::Float64, sUp::Float64, sDown::Float64, E_s::Float64)
    update = Vector{Float64}(undef,2)

    # Input voltage is from neuron 1 + some synapse input
    update[1] = ((E_L - v1) + (R_mI_e))    +    (R_mg_s * sDown * (E_s - v1))
    update[1] = update[1] * tau_m2

    # Input voltage is from neuron 2 + some synapse input
    update[2] = ((E_L - v2) + (R_mI_e))    +    (R_mg_s * sUp *   (E_s - v2))
    update[2] = update[2] * tau_m2

    update[1] += v1
    update[2] += v2

    return update
end

function updateSynapses(sUp::Float64, sDown::Float64)
    update = Vector{Float64}(undef,2)

    #Tau_s ds/dt = -s
    update[1] = -sUp
    update[1] = update[1] * tau_s

    update[2] = -sDown
    update[2] = update[2] * tau_s

    update[1] += sUp
    update[2] += sDown

    return update
end


function plotVoltages2(v1::Vector{Float64}, v2::Vector{Float64}, file::String)
    fig, ax = plt.subplots()
    ax.plot(v1, "b-")
    ax.plot(v2, "r--")
    plt.title("Modelled Pair of Linked Leaky Integrate and Fire Neurons\n")
    plt.xlabel("time/ms")
    plt.ylabel("Voltage/V")
    savefig(file)
    println("File saved to ", file, "\n")
end








function q2()
    voltage1 = zeros(1001)
    voltage1[1] = rand(Uniform(V_r2,V_t2))
    voltage2 = zeros(1001)
    voltage2[1] = rand(Uniform(V_r2,V_t2))

    synapseUp = 0.0::Float64
    synapseDown = 0.0::Float64

    #Storage vectors for return values from update functions
    update = Vector{Float64}(undef,2)
    updateS = Vector{Float64}(undef,2)

    for t in 2:1001
        update = updateVoltages2(voltage1[t-1], voltage2[t-1], synapseUp, synapseDown, 0.0)
        voltage1[t] = update[1]
        voltage2[t] = update[2]

        updateS = updateSynapses(synapseUp,synapseDown)
        synapseUp = updateS[1]
        synapseDown = updateS[2]

        if checkReset2(voltage1[t])
            voltage1[t] = V_r2
            synapseUp += P
        end
        if checkReset2(voltage2[t])
            voltage2[t] = V_r2
            synapseDown += P
        end
    end
    plotVoltages2(voltage1,voltage2, "pairEx.jpg")

    synapseUp = 0.0
    synapseDown = 0.0
    voltage1[1] = rand(Uniform(V_r2,V_t2))
    voltage2[1] = rand(Uniform(V_r2,V_t2))


    for t in 2:1001
        update = updateVoltages2(voltage1[t-1], voltage2[t-1], synapseUp, synapseDown, -0.08)
        voltage1[t] = update[1]
        voltage2[t] = update[2]

        if checkReset2(voltage1[t])
            voltage1[t] = V_r2
            synapseUp += P
        end
        if checkReset2(voltage2[t])
            voltage2[t] = V_r2
            synapseDown += P
        end
    end

    plotVoltages2(voltage1, voltage2, "pairIn.jpg")
end

function variableIUpdateVoltage(voltage::Float64, input_i::Float64)
    update = ((E_L - voltage) + (R_m*input_i))
    update = update * tau_m

    newVoltage = voltage + update

    return newVoltage
end


const global nA = 0.000000001

function plotSpikeCountsCurrent(spikes::Vector{Int64}, currents::Vector{Float64}, path::String)
    fig, ax = plt.subplots()
    ax.plot(currents,spikes, "b-")
    plt.title("Firing Rate of LINF neuron with variable Input Current\n")
    plt.xlabel("Input Current/A")
    plt.ylabel("Firing Rate/Hz")
    savefig(path)
    println("File saved to ", path, "\n")
end


function q3_3()
    voltage = zeros(1001)
    voltage[1] = V_r
    currentTrial = 1::Int64
    spikeCounts = zeros(Int64,31)
    currents = zeros(31)
    # 2:1001 to account for t[0] stored at voltage[1] due to Julia
    for i in 20:50
        currents[currentTrial] = i * 0.1 * nA
        for t in 2:1001
            voltage[t] = variableIUpdateVoltage(voltage[t-1], currents[currentTrial])
            if checkReset(voltage[t])
                voltage[t] = V_r
                spikeCounts[currentTrial] += 1
            end
        end
        currentTrial += 1
    end
    plotSpikeCountsCurrent(spikeCounts,currents,"firing.jpg")
end



# Call question functions
#q1()
#q2()
q3_3()
