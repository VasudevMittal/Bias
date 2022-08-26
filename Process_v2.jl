#Similar to Process.jl of V2. Differences are commented.
include("Simulation_v2.jl")

#Process function is done in Galactic Coordinates.
function process(i,t,nside)
    a = dataset(i)
    ra,dec = [a[i][1] for i in 1:length(a)],[a[i][2] for i in 1:length(a)]
    pu = scattomap(n[1:length(ra)],ra, dec,nside)
    b = doppler2.(pu,n2[1:length(pix)],pix)
    res = fit_dipole(b,nside,30)
    c = [res[2]/res[1],res[3]/res[1],res[4]/res[1]]
    d,e,f = bias(t,c),tht(t,c),norm(c)
    list = [res[2]/res[1],res[3]/res[1],res[4]/res[1],d,e,f]
    return list
end

i = 3200000
nside = 64
dx,dy,dz,b,o,m,filelist = [],[],[],[],[],[],[]
while length(m)<10000
    list = process(i,t,nside)
    push!(filelist,list),push!(dx,list[1]),push!(dy,list[2]),push!(dz,list[3]),push!(b,list[4]),push!(o,list[5]),push!(m,list[6])
end

coord = []
l = c2s.(dx,dy,dz)
lon,lat = [l[i][1] for i in 1:length(l)],[l[i][2] for i in 1:length(l)]
for i in 1:length(lon)
    if lon[i]>180
        lon[i]-=360
    end
    push!(coord,[lon[i],lat[i]])
end
writedlm("o.dat",coord,"\t")
writedlm("output.dat",filelist,"\t")