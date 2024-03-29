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
using Printf

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

function hopfield(a::Array{Int64}, b::Array{Int64}, c::Array{Int64}, w::Array{Float64})
    for i in 1:11
        for j in 1:11
            if i != j
                sum = 0

                sum += a[i] * a[j]
                sum += b[i] * b[j]
                sum += c[i] * c[j]

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
        sum = 0.0
        for j in 1:11
            if i != j
                sum += (w[i,j] * temp[j]) - theta
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

#Lines involving c were used to test cycling in the system on random test input
function evolve(pattern::Array{Int64}, w::Array{Float64}, f::IOStream)
    change = true
    #c = 0
    while change
        change = mp_update(pattern,w)
        seven_segment(pattern)
        seven_segment(f,pattern)
        qquad(f)
        print_number(f,energy(pattern,w))
        qquad(f)
        print_energy(energy(pattern,w))
        #c+=1
        #println(c)
        #if c > 15
        #    change = false
        #end
    end
end

#Add a function to calculate the energy of a configuration
# E = -1/2 SIGMA_ij x_i w_ij x_j

function energy(pattern::Array{Int64}, w::Array{Float64})
    sum = 0.0
    for i in 1:11
        for j in 1:11
            if i != j
                sum += pattern[i] * pattern[j] * w[i,j]
            end
        end
    end
    en = sum * (-1/2)
    return en
end

function print_energy(num)
    @printf "Energy of pattern : %f\n" num
end


f=open("./david_sharp.tex","w")
header(f,"David Sharp")

w = zeros(Float64, (11,11))

six=Int64[1,1,-1,1,1,1,1,-1,1,1,-1]
three=Int64[1,-1,1,1,-1,1,1,1,1,-1,-1]
one=Int64[-1,-1,1,-1,-1,1,-1,1,-1,-1,-1]

hopfield(six,three,one,w)

section(f,"Weight matrix")
matrix_print(f,"W",w)

seven_segment(three)
print_energy(energy(three,w))

seven_segment(six)
print_energy(energy(six,w))

seven_segment(one)
print_energy(energy(one,w))

cr(f)

#------------------
println("test1")
section(f,"Test Pattern 1")
test=Int64[1,-1,1,1,-1,1,1,-1,-1,-1,-1]

seven_segment(test)
seven_segment(f,test)

qquad(f)
print_number(f,energy(test,w))
qquad(f)
print_energy(energy(test,w))
#here the network should run printing at each step
evolve(test,w,f)
cr(f)

println("test2")
section(f,"Test Pattern 2")
test=Int64[1,1,1,1,1,1,1,-1,-1,-1,-1]

seven_segment(test)
seven_segment(f,test)

qquad(f)
print_number(f,energy(test,w))
qquad(f)
print_energy(energy(test,w))

#here the network should run printing at each step
evolve(test,w,f)
cr(f)

bottomer(f)
close(f)

#Random test pattern generation used to test the amount of cycles that exist
#in the system with three attractor states, found two pairs of states that
#would flip flop

#println("test3")
#test = rand([1,-1],11)
#seven_segment(test)
#print_energy(energy(test,w))
#evolve(test,w)
