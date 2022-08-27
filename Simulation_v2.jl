#Most of the functions are same as Simulation.jl However, boosting function is different. So commented functions reflect the changes.
using Healpix, AstroLib, StatsBase, Distributions, DelimitedFiles, CoordinateTransformations,StaticArrays,Rotations, LinearAlgebra, Statistics, Plots

function dataset(a)
    b = []
    c = [deg2rad(x) for x in vec(rand(Uniform(0,360),1,a))]
    mul = []
    while length(mul)<a
        push!(mul,rand([-1,1]))
    end
    d1 = [x for x in vec(rand(Uniform(0,1),1,a))].*mul
    d = asin.(d1)
    while length(b)<a
        push!(b,[c[rand(1:end)],d[rand(1:end)]])
    end
    return b
end

function s2c(alpha, delta)
    alpha,delta = deg2rad(alpha),deg2rad(delta)
    x = cos(alpha) * cos(delta)
    y = cos(delta) * sin(alpha)
    z = sin(delta)
    return [x, y, z]
end

function c2s(x, y, z)
    h = sqrt((x^2)+(y^2))
    alpha,delta = rad2deg(atan(y, x)),rad2deg(atan(z, h))
    if alpha<0;alpha+=360;end
    return [alpha, delta]
end

function scattomap(m, ra, dec , nside)
    data = ang2pixRing.(m,lat2colat.(dec), ra)
    hmap = counts(data,1:Healpix.nside2npix(nside))
    if length(hmap)==(Healpix.nside2npix(nside)+1)
        hmap=hmap[2:end]
    end
    return hmap
end

#Apart from the usual n containing the Healpix Resolution, n2 contains the nside for doppler2 function.
nside = 64
n,n2 = [],[]
while length(n)<10000000
    push!(n,Healpix.Resolution(nside)),push!(n2,nside)
end

#pix array contains the list of pixels. Useful for doppler2 function.
pix = [x for x in 1:12*nside*nside]

function pix2vec(ipix,nside)
    u = pix2angRing(Resolution(nside), ipix)
    mid = ang2vec(u[1],u[2])
    x,y,z = mid[1],mid[2],mid[3]
    return x,y,z
end   

#Pixelised boosting function.
function doppler2(p,q,r)#pixel,nside,ipix
    t = s2c(262.92103505720036,47.67440904728711)
    v,c,x,alpha = 369000,300000000,0.78,0.75
    beta = v/c
    cv = pix2vec(r,q)
    gamma = (1-(beta^2))^(-0.5)
    d = gamma*(1+(beta*cos(tht(cv,t))))
    e = d^(2+(x*(1+alpha)))
    f = e*p
    return f
end

function fit_dipole(m, nside, gal_cut=0,badval=UNSEEN)
    npix = size(m)[1]
    bunchsize = npix
    aa = zeros(4, 4)
    v = zeros(4)
    ipix = [x for x in 1:bunchsize]
    ipix = [ipix[i] for i in 1:length(ipix) if m[i]!=badval && isfinite(m[i])!=false]
    mid = pix2vec.(ipix,nside)
    x,y,z = [mid[i][1] for i in 1:length(mid)],[mid[i][2] for i in 1:length(mid)],[mid[i][3] for i in 1:length(mid)]
    if gal_cut > 0
        w = abs.(z) .>= sin(gal_cut * pi / 180)
        ipix = ipix[w]
        x = x[w]
        y = y[w]
        z = z[w]
    end
    mflatipix = [vec(m)[u] for u in ipix]
    aa[1, 1] += size(ipix)[1]
    aa[2, 1] += sum(x)
    aa[3, 1] += sum(y)
    aa[4, 1] += sum(z)
    aa[2, 2] += sum(x.^2)
    aa[3, 2] += sum(x.*y)
    aa[4, 2] += sum(x.*z)
    aa[3, 3] += sum(y.^2)
    aa[4, 3] += sum(y.*z)
    aa[4, 4] += sum(z.^2)
    v[1] += sum(mflatipix)
    v[2] += sum(mflatipix.*x)
    v[3] += sum(mflatipix.*y)
    v[4] += sum(mflatipix.*z)
    aa[1,2] = aa[2,1]
    aa[1,3] = aa[3,1]
    aa[1,4] = aa[4,1]
    aa[2,3] = aa[3,2]
    aa[2,4] = aa[3,2]
    aa[3,4] = aa[4,3]
    res = inv(aa)*v
    return res
end

v,c,x,alpha = 369000,300000000,0.78,0.75
r = (2+(x*(1+alpha)))*(v/c)
t = r*s2c(262.92103505720036,47.67440904728711)

function bias(t,b)
    return norm(t+b)/norm(t)
end

function tht(t,u)
    return acos(dot(t,u)/(norm(t)*norm(u)))
end

