# just for fun script to calculate the trajectorie of a car based on the state space description of the single track model

# dβ/dt = −(cv+ch)/(m*v)*β + (m*v^2−(ch*lh−cv*lv))/(m*v^2)*dψ/dt − cv/(m*v)*δH/iS
# d^2ψ/dt^2 = −(ch*lh−cv*lv)/θ*β − (ch*lh^2+cv*lv^2)/(θ·v)*dψ/dt + (cv*lv)/θ*δH/iS
# Source: https://de.wikipedia.org/wiki/Einspurmodell

using Plots

include("NumericalMethods.jl")


t_intervall = 0.01
t_start = 10
t_stop = 30
time::Vector{Float64} = t_start:t_intervall:t_stop
δH::Vector{Float64} = pi*time/10
v = 30/3.6


struct Car
    m::Int16   # mass
    θz::Float64  # moment of inertia about z-Axis

    lf::Float64  # distance front axle to center of gravity
    lr::Float64  # distance rear axle to center of gravity

    cf::Int32   # side slip stiffness front
    cr::Int32   # side slip stiffness rear
    iS::Int8  # steering ration

    #Matricies for state space description 
    a::Matrix{Float64}
    b::Vector{Float64}
    c::Matrix{Int8}
    d::Vector{Int8}

    function Car(m,θz,lf,lr,cf,cr,iS)
        #defining state space matricies for single track model
        a = [−(cf+cr)/(m*v) (m*v^2−(cr*lr−cf*lf))/(m*v^2);−(cr*lr−cf*lf)/θz −(cr*lr^2+cf*lf^2)/(θz*v)]

        b = [-cf/(m*v);(cf*lf)/θz]

        c = [1 0;0 1]

        d = [0;0]
        
        new(m,θz,lf,lr,cf,cr,iS, a, b, c, d)
    end

    
    
end


car = Car(1550, 2800, 1.344, 1.456, 75000, 150000, 16)


u::Vector{Float64} = δH/car.iS #input vector of state space description
x0::Vector{Float64} = [0;0] #initial condition [β dψ]

carModel(x,u) = car.a*x+car.b*u # state space equation to pass over to solver
solution::Matrix{Float64} = rungeKutta(carModel, x0, u, t_intervall) # [β1 ... βn;dψ1 ... dψn]

solution = solution[:,1:end-1] #cutting solution to size


# extracting values from solution
β::Vector{Float64} = solution[1,:] # [β1;...;βn]
dψ::Vector{Float64} = solution[2,:] # same for dψ
ψ::Vector{Float64} = trapez(dψ, t_intervall) #integrate to get yaw angle

# position based on angle and current speed
x::Vector{Float64} = trapez(cos.(β + ψ)*v, t_intervall)
y::Vector{Float64} = trapez(sin.(β + ψ)*v, t_intervall)

# position of the front and rear axle
car_orientation_vector::Matrix{Float64} = [cos.(ψ) sin.(ψ)]
front_tire::Matrix{Float64} = [x y] + car.lf*car_orientation_vector # [x1 y1;...;xn yn]
rear_tire::Matrix{Float64} = [x y] - car.lr*car_orientation_vector






                        #= Visualization =#


#for correct temporal depiction, the framerate has to match the sample rate
animation_speed = 1
framerate = 30
frame_intervall::Int64 = round(animation_speed/framerate/(time[2]-time[1])) # round might 
frame_total = length(x)

#exaggerated tire size
tire_width::Float32 = 0.4
tire_diameter::Float32 = 0.7

trace_length = 1000000 #how many frames the trace should last

@userplot SingleTrackModel
@recipe function f(stm::SingleTrackModel)
    front::Matrix{Float32}, rear::Matrix{Float32} = stm.args
    legend --> false

    aspect_ratio --> 1
    label --> false
    #ylims --> (floor(minimum(y)-5),floor(maximum(y)+5))
    

    [front[1], rear[1]], [front[2], rear[2]]
end



@userplot Tire
@recipe function f(axle::Tire)
    coordinates, α  = axle.args
    linewidth --> 2
    linecolor --> :black

    # calculating the corners of the tire 
    tire_orientation_vector::Matrix{Float32} = [cos(α) sin(α)] # the following calculations result in 5x2 Matrix
    tire_bounds::Matrix{Float32} = ones(5)*coordinates + [1;1;-1;-1;1]*tire_orientation_vector*tire_diameter/2 + [1;-1;-1;1;1]*tire_orientation_vector*[0 1;-1 0]*tire_width/2    
    
    # the original less elegant calculation with much better readability
    # tire_front_middle =         coordinates + tire_orientation_vector*tire_diameter/2
    # tire_rear_middle =          coordinates - tire_orientation_vector*              tire_diameter/2

    # tire_front_left =           tire_front_middle + tire_orientation_vector*[0 1;-1 0]*tire_width/2
    # tire_front_right =          tire_front_middle - tire_orientation_vector*[0 1;-1 0]*tire_width/2
    # tire_rear_right =           tire_rear_middle - tire_orientation_vector*[0 1;-1 0]*tire_width/2
    # tire_rear_left =            tire_rear_middle + tire_orientation_vector*[0 1;-1 0]*tire_width/2
    # tire_front_left =           tire_front_middle + tire_orientation_vector*[0 1;-1 0]*tire_width/2

    # tire_bounds =               [tire_front_left;tire_front_right;tire_rear_right;tire_rear_left;tire_front_left]

    tire_bounds[:,1], tire_bounds[:,2]
end


anim = @animate for i ∈ 1:frame_intervall:frame_total
    #calls recipe function "SingleTrackModel", provides coordinates of both wheels
    front::Matrix{Float32}, rear::Matrix{Float32} = transpose(front_tire[i,:]), transpose(rear_tire[i,:])

    singletrackmodel(front, rear)# body of stm-depiction
    tire!(front, ψ[i] + δH[i]/car.iS) #front tire
    tire!(rear, ψ[i]) #rear tire

    #trajectory trace
    if trace_length==0
    elseif trace_length > i
        plot!(x[begin:i],y[begin:i])
        #plot!(front_tire[begin:i,1],front_tire[begin:i,2])
        #plot!(rear_tire[begin:i,1],rear_tire[begin:i,2])
    else
        plot!(x[i-trace_length:i],y[i-trace_length:i])
        #plot!(front_tire[i-trace_length:i,1],front_tire[i-trace_length:i,2])
        #plot!(rear_tire[i-trace_length:i,1],rear_tire[i-trace_length:i,2])
    end
end
gif(anim, "stm.gif", fps = framerate)