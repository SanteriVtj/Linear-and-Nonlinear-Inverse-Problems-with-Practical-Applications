using LinearAlgebra, Trapz, Plots, Interpolations, Distributions, Test, SparseArrays

function add_pct_noise(x,sigma)
    rand.(Normal.(x ,sigma.*abs.(uData)),1) |>
        Iterators.flatten |>
        collect
end

sigma = .01
L = π
D = 1
T = 1
xdim = 24
tdim = 108

x = LinRange(0, L, xdim)
t = LinRange(0, T, tdim)
dx = x[2]-x[1]
dt = t[2]-t[1]

@test dt/dx^2>1/2

# Initial heat distribution
Z = zeros(tdim,xdim)
Z[1,:] .= 10*sin.(2*x)   # u(x,0) = 10sin(2x)
Z[:,1] .= 0                 # u(0,t) = 0
Z[:,xdim] .= 0              # u(L,t) = 0

# Simulate the heat equation up to T
v = 2:xdim-1
for time in 2:tdim
    time_m1 = time-1
    Z[time,v] .= Z[time_m1,v] .+ dt/dx^2*(Z[time_m1,v.+1].-2*Z[time_m1,v].+Z[time_m1,v.-1])
end

ZM = zeros(xdim, tdim)
ZM[:,1] .= 10*sin.(2*x)   # u(x,0) = 10sin(2x)
ZM[1,:] .= 0                 # u(0,t) = 0
ZM[xdim,:] .= 0              # u(L,t) = 0

S = D*dt/dx^2
n = length(x)-2
A = (1-2*S)*Matrix(I(n))
A[diagind(A, 1)] .= S
A[diagind(A, -1)] .= S

for time in 2:tdim
    ZM[v, time] .= A*ZM[v,time-1]
end

uData = Z[tdim,:]
uDataM = ZM[:,tdim]

uDataNoisy = rand.(Normal.(uData ,sigma.*abs.(uData)),1) |>
    Iterators.flatten |>
    collect

# Exact solution
function u_exact(x;t=T)
    10*exp(-4*t)*sin(2*x)
end

sn = maximum(u_exact.(x))
supnorm = norm(Z[tdim,:].-u_exact.(x),Inf)

plot(x, Z[1,:])
plot!(x, Z[tdim,:])

plot(t,x,Z',st=:surface,camera=(-120,30))

# Inverse solver
function solve_inverse(u;exact=true,n_it=4)
    xn = x[2:xdim-1]
    u = u[2:end-1]
    ss = tdim-n_it+1
    
    A_inv = A\I
    back_step = A_inv*u
    p = plot()
    plot!(p, xn, back_step,label="1")
    for i in 1:tdim-ss
        back_step .= A_inv*back_step
        plot!(p,xn, back_step, label="$(i+1)")
    end
    if exact
        plot!(p,xn,u_exact.(xn,t=dt*ss),label="exact",color=:red)
    end
    p
end

solve_inverse(uData)

ZI = zeros(xdim, tdim)
ZI[:,1] .= 10 .*(x.>=π/4).*(x.<=3*π/4)   # u(x,0) = 10sin(2x)
ZI[1,:] .= 0                 # u(0,t) = 0
ZI[xdim,:] .= 0    
for time in 2:tdim
    ZI[v, time] .= A*ZI[v,time-1]
end
ud = ZI[:,tdim]
ud_noisy = add_pct_noise(ud,sigma)

solve_inverse(ud, exact=false)

solve_inverse(ud_noisy,exact=false,n_it=3)