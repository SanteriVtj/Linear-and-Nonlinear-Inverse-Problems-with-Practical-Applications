using LinearAlgebra, Trapz, Plots, Interpolations, Distributions, Test, SparseArrays, ImagePhantoms, RadonKA

n = 50
sl = shepp_logan(n,n)
xrange,yrange = LinRange(0,1,n),LinRange(0,1,n)

heatmap(xrange,yrange,sl',color=:greys)

# θ = collect(range(0f0, 2f0π, 180)[begin:end-1])
θ = collect(range(0f0, 2f0π, 180))

R = RadonKA.radon(sl, θ)
heatmap(R)
bp = RadonKA.backproject(R, θ)
heatmap(bp)