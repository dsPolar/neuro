I had some issues with load.jl and poisson.jl, generally just using functions
from modules without stating that in code. I also had some issues with
the skeleton using Vector{Float64}(0) being used to try and initialise a 0-size
vector so altered that.

As a result I have also attached the modified load.jl and poisson.jl files
so that my code will run if any problems arise.
