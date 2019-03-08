#Create a weight matrix of w_ij
#Fix the weight values to store the three patterns
#Write a function to evolve the network according to MP formula
#Done synchronously so every node is updated at each time step
#For each of the two test patterns, evolve the patterns until they
#stop changing, printing the pattern at each step

#Masters step
#Add a function to calculate the energy of a configuration
# E = -1/2 SIGMA_ij x_i w_ij x_j
#Print out the energy of teh three learned patterns
#Of the the tests patterns
#Any patterns formed as the patterns are updated

include("submission.jl")

const global theta = 0

#Takes an 11 component vector and prints the corresponding seven segment
function seven_segment(pattern::Array{Int64})

    function to_bool(a::Int64)
        if a==1
            return true
        end
        false
    end

    function hor(d)
        if d
            println(" _ ")
        else
            println("   ")
        end
    end


    function vert(d1,d2,d3)
        word=""::String
        if d1
            word="|"
        else
            word=" "
        end
        if d3
            word=string(word,"_")
        else
            word=string(word," ")
        end
        if d2
            word=string(word,"|")
        else
            word=string(word," ")
        end
        println(word)

    end

    pattern_b=map(to_bool,pattern)

    hor(pattern_b[1])
    vert(pattern_b[2],pattern_b[3],pattern_b[4])
    vert(pattern_b[5],pattern_b[6],pattern_b[7])

    number=0
    for i in 1:4
        if pattern_b[7+i]
            number+=2^(i-1)
        end
    end

    println(number)

end

#Patterns is a single array containing all three attractor states
#Accessed here by using state offset + 11 * state number
function hopfield(patterns::Array{Int64}, w::Array{Float64})
    for i in 1:11
        for j in 1:11
            if i != j
                sum = 0
                for p in 0:2
                    sum += patterns[p*11+i] * patterns[p*11+j]
                end
                w[i,j] = (1/3) * sum
            end
        end
    end
end

function mp_update(pattern::Array{Int64}, w::Array{Float64})
    function g(z)
        val = -1
        if z > 0
            val = 1
        end
        return val
    end

    change = false
    sum = 0.0
    temp = copy(pattern)

    for i in 1:11
        for j in 1:11
            if i != j
                sum += w[i,j] * temp[j] - theta
            end
        end
        pattern[i] = g(sum)
        if !change
            if(pattern[i] != temp[i])
                change = true
            end
        end
    end
    return change
end

function evolve(pattern::Array{Int64}, w::Array{Float64})
    change = true

    while change
        change = mp_update(pattern,w)
        seven_segment(pattern)
    end
end



six=Int64[1,1,-1,1,1,1,1,-1,1,1,-1]
three=Int64[1,-1,1,1,-1,1,1,1,1,-1,-1]
one=Int64[-1,-1,1,-1,-1,1,-1,1,-1,-1,-1]

seven_segment(three)
seven_segment(six)
seven_segment(one)

#------------------
w = zeros(Float64, (11,11))

hopfield([six;three;one],w)

println("test1")

test=Int64[1,-1,1,1,-1,1,1,-1,-1,-1,-1]

seven_segment(test)

#here the network should run printing at each step
evolve(test,w)

println("test2")

test=Int64[1,1,1,1,1,1,1,-1,-1,-1,-1]

seven_segment(test)

#here the network should run printing at each step
evolve(test,w)
