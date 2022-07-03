using Healpix, AstroLib, StatsBase, Distributions, DelimitedFiles, CoordinateTransformations,StaticArrays,Rotations, LinearAlgebra

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

function rotor(d,a,b,c)
    x = s2c(d[1],d[2])
    ma = transpose(RotZXZ(deg2rad(90+a),deg2rad(90-b),deg2rad(-(c-90))))
    y = ma*x
    return c2s(y[1],y[2],y[3])
end


function rotorback(d,a,b,c)
    x = s2c(d[1],d[2])
    mb = transpose(RotZXZ(deg2rad(c-90),deg2rad(b-90),deg2rad(-90-a)))
    y = mb*x
    return c2s(y[1],y[2],y[3])
end

function doppler(d)#radec#degrees#v#m/s#arrayradec#radians
    a,b,p,c,v = 167,-7,5,3*(10^8), 370000
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

function galcut(a,b,c)
    if abs(euler(rad2deg(a),rad2deg(b),1)[2])>c
        return [a,b]
    end
end

function sgalcut(a,b,c)
    u1,v1 = euler(rad2deg(a),rad2deg(b),1)[1],euler(rad2deg(a),rad2deg(b),1)[2]#rd
    if abs(gal2sgal(u1,v1))>c
        return [a,b]
    end
end

function scattomap(m, ra, dec , nside)
    data = ang2pixRing.(m,lat2colat.(dec), ra)
    hmap = counts(data,1:Healpix.nside2npix(nside))
    if length(hmap)==(Healpix.nside2npix(nside)+1)
        hmap=hmap[2:end]
    end
    return hmap
end

nside = 64
n = []
while length(n)<10000000
    push!(n,Healpix.Resolution(nside))
end