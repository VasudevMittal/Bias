using Healpix, AstroLib, StatsBase, Distributions, DelimitedFiles, CoordinateTransformations,StaticArrays,Rotations, LinearAlgebra, Statistics, Plots

#Function to convert from Galactic to Super-Galactic Coordinates
R=[-0.735743    0.677261   0.0
   -0.0745538  -0.0809915  0.993923
    0.673145    0.731271   0.110081]
function gal2sgal(l,b)#dd
    SGB=float(b)
    l=deg2rad(l);b=deg2rad(b);
    X_sgal=R*[sin(pi/2-b)*cos(l),sin(pi/2-b)*sin(l),cos(pi/2-b)]
    SGB=rad2deg(atan(X_sgal[3],sqrt(X_sgal[1]^2+X_sgal[2]^2)))    
    return SGB     
end

#Function to generate a full sky map of "a" sources.
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

#Function to convert from lat-lon to cartesian coordinates.
function s2c(alpha, delta)
    alpha,delta = deg2rad(alpha),deg2rad(delta)
    x = cos(alpha) * cos(delta)
    y = cos(delta) * sin(alpha)
    z = sin(delta)
    return [x, y, z]
end

#Function to convert from cartesian coordinates to lat-lon.
function c2s(x, y, z)
    h = sqrt((x^2)+(y^2))
    alpha,delta = rad2deg(atan(y, x)),rad2deg(atan(z, h))
    if alpha<0;alpha+=360;end
    return [alpha, delta]
end

#Function to rotate the sky from celestial to native coordinates using the ZXZ euler rotation with angles of a,b and c.
function rotor(d,a,b,c)
    x = s2c(d[1],d[2])
    ma = transpose(RotZXZ(deg2rad(90+a),deg2rad(90-b),deg2rad(-(c-90))))
    y = ma*x
    return c2s(y[1],y[2],y[3])
end

#Function for inverse rotation of sky from native to celestial coordinates.
function rotorback(d,a,b,c)
    x = s2c(d[1],d[2])
    mb = transpose(RotZXZ(deg2rad(c-90),deg2rad(b-90),deg2rad(-90-a)))
    y = mb*x
    return c2s(y[1],y[2],y[3])
end

#Function for giving a doppler boost to the sources in the direction of CMB dipole. Here, the dipole is located at (167,-7) in RA-DEC and the velocity of the observer is 369km/s. 5 is the pole_longitude.
function doppler(d)
    a,b,p,c,v = 167,-7,5,3*(10^8), 369000
    gamma = 1/sqrt(1-((v/c)^2))
    ra,dec = deg2rad(a),deg2rad(b)
    dangle = acos(cos(d[2])*cos(dec)*cos(d[1] - ra)+sin(d[2])*sin(dec))
    tanphi = sin(dangle)/(gamma*(cos(dangle)-v/c))
    lon, lat = rotor([rad2deg(d[1]), rad2deg(d[2])],a,b,p)
    phi = rad2deg(atan(tanphi))
    phi = phi+(-1*sign(phi) + 1)*90
    diff = rad2deg(dangle) - phi
    latab = lat - diff
    k,l = rotorback([lon, latab],a,b,p)
    return [deg2rad(k),deg2rad(l)]
end

#Function to apply a galactic cut to a sky in equatorial coordinates within galactic latitudes of +-c degrees about galactic equator.
function galcut(a,b,c)
    if abs(euler(rad2deg(a),rad2deg(b),1)[2])>c
        return [a,b]
    end
end

#Function to apply a super-galactic cut to a sky in equatorial coordinates within galactic latitudes of +-c degrees about super-galactic equator.
function sgalcut(a,b,c)
    u1,v1 = euler(rad2deg(a),rad2deg(b),1)[1],euler(rad2deg(a),rad2deg(b),1)[2]#rd
    if abs(gal2sgal(u1,v1))>c
        return [a,b]
    end
end

#Function to convert the sky map into a Healpix map on ring order with a given nside.
function scattomap(m, ra, dec , nside)
    data = ang2pixRing.(m,lat2colat.(dec), ra)
    hmap = counts(data,1:Healpix.nside2npix(nside))
    if length(hmap)==(Healpix.nside2npix(nside)+1)
        hmap=hmap[2:end]
    end
    return hmap
end

#Creation of this array is useful for further processing.
nside = 64
n = []
while length(n)<10000000
    push!(n,Healpix.Resolution(nside))
end

#Function to convert a Healpix Pixel of Ring Order to a vector in cartesian coordinates
function pix2vec(ipix,nside)
    u = pix2angRing(Resolution(nside), ipix)
    mid = ang2vec(u[1],u[2])
    x,y,z = mid[1],mid[2],mid[3]
    return x,y,z
end   

#Function to implement fit_dipole routine of healpy. Presently, if a mask is to be applied, the pixel values should be 0. Maps with bad pixels or nside>>128 are not supported.
function fit_dipole(m, nside, gal_cut=0)
    npix = size(m)[1]
    bunchsize = npix
    aa = zeros(4, 4)
    v = zeros(4)
    ibunch = 0
    ipix = [x for x in 1:bunchsize]
    ipix = ipix[1:maximum(ipix)]
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

#Function to calculate the bias.
function bias(t,b)
    return norm(t+b)/norm(t)
end

#Function to calculate the offset from a given dipole.
function tht(t,u)
    return acos(dot(t,u)/(norm(t)*norm(u)))
end

