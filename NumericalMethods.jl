# I might reuse these, so I defined them in a seperate file


function rungeKutta(f,x0,u,h)
    #Solution to differential equation f(x,u) using Runge Kutta

    #setup of result vektor
    x = zeros(Float64, (length(x0), length(u)+1))
    #setting initial condition
    x[:,1]=x0

    #in this case:
    #x[:,i] == yi
    #x[:,i+1] == yi+1

    for i in eachindex(u)
        #calculating intermediate steps
        k1 = f(x[:,i],u[i])
        k2 = f((x[:,i]+h*k1/2),u[i])
        k3 = f((x[:,i]+h*k2/2),u[i])
        k4 = f((x[:,i]+h*k3),u[i])

        #calculating yk+1
        x[:,(i+1)] = x[:,i]+h*(k1+k2+k3+k4)/6
    end
    return x 
end

function trapez(y::Vector{Float64},x::Vector{Float64})
    #integrate y over x using trapezoidal approximation
    n=length(y)
    ∫ydx = zeros(Float64, n)
    h = x[2]-x[1] # assuming equidistant intervalls

    for i in 1:length(∫ydx)-1
        ∫ydx[i+1] = (y[i] + y[i+1])*h/2 + ∫ydx[i]
    end

    return ∫ydx
end


function trapez(y::Vector{Float64},h::Float64)
    #integrate using trapezoidal approximation, this time the intervall of x is given as an argument
    n=length(y)
    ∫ydx = zeros(Float64, n)

    for i in 1:length(∫ydx)-1
        ∫ydx[i+1] = (y[i] + y[i+1])*h/2 + ∫ydx[i]
    end

    return ∫ydx
end

#using Plots
#x::Vector{Float64} = 0:10
#y = ones(Float64, 11)
#Y = trapez(y, x)
#plot(x,Y)