include("Simulation.jl")

function process(i,t,nside)
    a = dataset(i)
    dop = doppler.(a)
    ra,dec = [dop[i][1] for i in 1:length(dop)],[dop[i][2] for i in 1:length(dop)]
    gal= galcut.(ra,dec,30)
    gal=gal[gal.!=nothing]
    galra,galdec = [gal[i][1] for i in 1:length(gal)],[gal[i][2] for i in 1:length(gal)]
    b = scattomap(n[1:length(galra)],galra, galdec,nside)
    res = fit_dipole(b,nside)
    c = [res[2]/res[1],res[3]/res[1],res[4]/res[1]]
    d,e,f = bias(t,c),tht(t,c),norm(c)
    list = [res[2]/res[1],res[3]/res[1],res[4]/res[1],d,e,f]
    return list
end

i = 3200000
r,a,d = 0.007,167,-7
t = r*s2c(a,d)
nside = 64
dx,dy,dz,b,o,m,filelist = [],[],[],[],[],[],[]
j = 0
while j<10000
    list = process(i,t,nside)
    push!(filelist,list),push!(dx,list[1]),push!(dy,list[2]),push!(dz,list[3]),push!(b,list[4]),push!(o,list[5]),push!(m,list[6])
    j+=1
end

m1,o1 = [0.007],[0]
plot(m,o,seriestype = :scatter,xaxis=:log,grid=true,minorgrid=true,label="Random Dipoles")
plot!(m1,o1,seriestype = :scatter,xaxis=:log,grid=true,minorgrid=true,label="True Dipoles")
xlabel!("Magnitude")
ylabel!("Offset")
title!("Magnitude vs Offset (30 deg Galactic Cut)")
savefig("ScatterPlot.png")

coord = []
l = c2s.(dx,dy,dz)
lon,lat = [l[i][1] for i in 1:length(l)],[l[i][2] for i in 1:length(l)]
gal = euler.(lon,lat,1)
galra,galdec = [gal[i][1]+180 for i in 1:length(gal)],[gal[i][2] for i in 1:length(gal)]
for i in 1:length(galra)
    if galra[i]>180
        galra[i]-=360
    end
    push!(coord,[galra[i],galdec[i]])
end
writedlm("o.dat",coord,"\t")
writedlm("output.dat",filelist,"\t")