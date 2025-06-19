using LinearAlgebra, Trapz, Plots, Interpolations, Distributions

# Define the point spread function
function PSF(x,a)
    if ((a<=0) | (a>=1/2))
        error("Parameter a should satisfy 0<a<1/2")
    end
    (x.+a).^2 .*(x.-a).^2
end

function target(x)
    f = zeros(size(x))

    x = x.-floor.(x)

    f[((x .>= 0.12) .& (x .<= 0.15))] .= 1.5
    f[((x .>= 0.2) .& (x .<= 0.25))] .= 1.3
    f[((x .>= 0.75) .& (x .<= 0.85))] .= 1
    ind = (x .>= 0.35) .& (x .<= 0.55)
    f[ind] .= 5*(x[ind].-0.35)
    return f
end


# Problem generation
a = .01

# Points for function evaluation
Nxx = 2000
xx = LinRange(0,1,Nxx)

# Points for integration
Nxxp = 1000
xxp = LinRange(-a,a,Nxxp)
Dxxp = xxp[2]-xxp[1]

# Evaluate point spread function
psf = zeros(Nxxp)
ind = abs.(xxp).<a
psf[ind] .= PSF(xxp[ind], a)
Ca = trapz(1:length(psf), psf)
Ca = 1/(Dxxp*Ca)
psf = Ca*psf

result = zeros(Nxx)
# Compute the convolution (Def. 2.1.1)
for i in 1:Nxx
    # ψ(x-x')
    targ = target(xx[i].-xxp)
    # ∫f(x')ψ(x-x')dx'
    result[i] = Dxxp*trapz(1:length(targ),psf.*targ)
end

plot(xx, result, label="covoluted")
plot!(xx, target(xx), label="true")

# Discretizised Problem
n = 64
x = collect(0:n-1)/n
Dx = x[2]-x[1]

sigma = 0.05

# Interpolate values
int = Interpolations.scale(interpolate(result, BSpline(Cubic(Line(OnGrid())))), xx)
data = int(x)

# Add noise
noisy_data = rand.(Normal.(data ,sigma.*abs.(data)),1) |>
    Iterators.flatten |>
    collect


plot(x, data)
scatter!(x, noisy_data)

max_error = maximum(abs.(data.-noisy_data)/maximum(abs.(data)))
println("Maximum error $max_error")

# Inverse crimes
nPSF = ceil(a/Dx)
xPSF = (-nPSF:nPSF)*Dx

PSF_data = zeros(size(xPSF))
ind = abs.(xPSF).<a
PSF_data[ind] .= PSF(PSF_data[ind],a)
Ca = 1/(Dx*trapz(1:length(PSF_data), PSF_data))
PSF_data = Ca*PSF_data

function convmtx(h::Vector, n::Int)
    k = length(h)
    m = n + k - 1
    H = zeros(eltype(h), m, n)

    for col in 1:n
        for row in 1:k
            if col + row - 1 <= m
                H[col + row - 1, col] = h[row]
            end
        end
    end
    H = H'
    # Correct for boundary condition
    m = floor(Int((length(h)-1)/2))
    H2 = H[:,m+1:end-m]
    for _ in 1:m
        H2[1:m,end-m+1:end] = H[1:m,1:m]
        H2[end-m+1:end,1:m] = H[end-m+1:end,end-m+1:end]
    end

    return H2
end

A = Dx*convmtx(PSF_data, n)

# Inverse crimed data
f = target(x)
mIC = A*f

plot(x, mIC)
scatter!(x, mIC)
scatter!(x, noisy_data)
plot!(x, target(x))

# A^-1*(A*f)
plot(x,(A\I)*mIC)
# A^-1*̂Af
plot!(x,(A\I)*data)
# A^-1*̂A(f+ϵ)
plot!(x,(A\I)*noisy_data)
