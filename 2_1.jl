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
Ca = trapz(xxp, psf)
Ca = 1/(Dxxp*Ca)
psf = Ca*psf

result = zeros(Nxx)
# Compute the convolution (Def. 2.1.1)
for i in 1:Nxx
    # ψ(x-x')
    targ = target(xx[i].-xxp)
    # ∫f(x')ψ(x-x')dx'
    result[i] = Dxxp*trapz(xxp,psf.*targ)
end

plot(xx, result, label="covoluted")
plot!(xx, target(xx), label="true")

# Discretizised Problem
n = 64
x = collect(0:n-1)/n
Dx = x[2]-x[1]

noise_level = 0.05

# Interpolate values
int = scale(interpolate(result, BSpline(Cubic(Line(OnGrid())))), xx)

# Add noise
noisyx = 

plot(x, int(x))
