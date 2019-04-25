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
using Statistics
using PyPlot

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

end

#Check if voltage is above V_t
function checkReset(voltage::Float64)
    if voltage > V_t
        return true
    else
        return false
    end
end

#Return reset voltage
function resetVoltage()
    return V_r
end

function plotVoltage(voltages::Vector{Float64})
    plot(voltages)
    plt.title("Modelled Leaky Integrate and Fire Neuron\n")
    plt.xlabel("Time step/ms")
    plt.ylabel("Membrane Voltage/V")
    savefig("linf.jpg")
    println("File saved to linf.jpg")
end

function q1()
    voltage = zeros(1001)
    for t in 1:1000
        voltage[t] = voltageUpdate(voltage[t-1])
        if checkReset(voltage[t])
            voltage[t] = resetVoltage()
        end
    end

end
