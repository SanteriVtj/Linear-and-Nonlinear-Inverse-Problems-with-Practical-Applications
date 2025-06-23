using LinearAlgebra, Trapz, Plots, Interpolations, Distributions, Test, SparseArrays, ImagePhantoms, RadonKA

function dist(x::Base.AbstractCartesianIndex, y::Base.AbstractCartesianIndex)
    return norm(x.I.-y.I)
end

function three_circles(;n=50, centers=[CartesianIndex(25,10), CartesianIndex(25,25), CartesianIndex(25,40)], R=[7,3,7],values=[.6,.7,.8])
    @assert length(centers) == length(R)
    X = zeros(n,n)
    indices = findall(x->true, X)
    for (rind,i) in enumerate(centers)
        X[dist.(i,indices).<=R[rind]] .= values[rind]
    end
    return X
end


n = 50
sl = shepp_logan(n,n)
xrange,yrange = LinRange(0,1,n),LinRange(0,1,n)

heatmap(xrange,yrange,sl',color=:greys)

θ = collect(range(0f0, 2f0π, 360))

R = RadonKA.radon(sl, θ)
heatmap(R)
bp = RadonKA.backproject(R, θ)
heatmap(bp)


tc = three_circles()
heatmap(tc)

t = [0, π/4, π/2, 3π/2, π]
tcRadon = RadonKA.radon(tc, θ)
heatmap(tcRadon)
tcBackRadon = RadonKA.backproject(tcRadon, θ)
heatmap(tcBackRadon)