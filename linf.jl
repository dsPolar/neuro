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

#Return reset voltage
function resetVoltage()
    return V_r
end


#Use pyplot to plot voltages vector and save to file
function plotVoltage(voltages::Vector{Float64})
    plot(voltages)
    plt.title("Modelled Leaky Integrate and Fire Neuron\n")
    plt.xlabel("time/ms")
    plt.ylabel("Voltage/V")
    savefig("linf.jpg")
    println("File saved to linf.jpg")
end



#Model a single linf neuron and plot voltage against time
function q1()
    voltage = zeros(1001)
    voltage[1] = V_r
    # 2:1001 to account for t[0] stored at voltage[1] due to Julia
    for t in 2:1001
        if checkReset(voltage[t-1])
            voltage[t] = voltageUpdate(resetVoltage())
        else
            voltage[t] = voltageUpdate(voltage[t-1])
        end

    end
    plotVoltage(voltage)
end





# Call question functions
q1()
