### Script for Running GrowClust: Julia

### Import packages

# Packages
using Distributed
using Printf
using DataFrames
using CSVFiles
using StatsKit
using Random
using Dates
using Interpolations
using SharedArrays

# GrowClust1D
using GrowClust1D

### Define Algorithm Parameters

# ------- GrowClust algorithm parameters -------------------------------
const distmax = 5.0          # maximum catalog(input) distance to join clusters (km)
const distmax2 = 3.0         # maximum relocated distance to join clusters (km)
const hshiftmax = 2.0        # maximum permitted horizontal cluster shifts (km)
const vshiftmax = 2.0        # maximum permitted vertical cluster shifts (km)
const rmedmax = 0.05         # maximum median absolute tdif residual to join clusters
const maxlink = 10           # use 10 best event pairs to relocate (optimize later...)
const nupdate = 10000         # update progress every nupdate pairs - NEW
   
# ------- Relative Relocation subroutine parameters -------------
const boxwid = 3. # initial "shrinking-box" width (km)
const nit = 15 # number of iterations
const irelonorm = 1 # relocation norm (L1 norm=1, L2 norm=2, 3=robust L2)
const tdifmax = 30. # maximum differential time value allowed (for error-checking on xcor data input)
const torgdifmax = 10.0 # maximum origin time adjustment (for robustness)
   
# -------- Bootstrap resampling parameters -------------------
const iseed = 0 # random number seed
#  1: Resample each event pair independently (each pair always has the same # of picks in each resample)'
#  2: Resample the entire data vectors at once (event pairs may have different # of picks in each resample)
        
# ------- Velocity model parameters (added 04/2018) -------------------------
const vzmodel_type = 1 # velocity model type: 1 = flat earth, (Z,Vp,Vs)
                 #                  or  2 = radial, (R,Vp,Vs): 
                 #                note: option 2 has not been extensively tested
const shallowmode = "flat" # option for how to treat shallow seismicity
                       # flat treats negative depths as zero depth
                       # reflect treats negative depths as equivalent to -depth; ok for true elevations

# ------- Geodetic parameters -------------
const degkm = 111.1949266


### Read Input File

# read input file
println("\nReading input file: ",ARGS[1])
infile_ctl = ARGS[1]
inpD = read_gcinp(infile_ctl)

# output paths
println("Assigning output directories:")
for fkey in collect(keys(inpD))
    if startswith(fkey,"fout")
        outdir = join(split(inpD[fkey],"/")[1:end-1],"/")
        println(fkey, " => ",outdir)
        mkpath(outdir)
    end
end

### Update fields

# set default Vp/Vs
if inpD["vpvs_factor"] < 0.01
    inpD["vpvs_factor"] = sqrt(3.0)
end

# define minimum ray parameter
if inpD["rayparam_min"] < 0.0
    inpD["plongcutP"] = 1.0/7.5 # assumes a sub-moho Vp of 7.5 km/s
    inpD["plongcutS"] = sqrt(3)/7.5 # Vs = Vp/sqrt(3) by default
else
    inpD["plongcutP"] = rayparam_min # user defined value
    inpD["plongcutS"] = rayparam_min*vpvs_factor # user defined value
end

# interpolation spacing
inpD["itp_dz"] = inpD["tt_ddep"]



### Check Input parameters
input_ok = check_gcinp(inpD)
if input_ok
    println("Input parameters are ok!")
else
    println("ERROR: FIX INPUT PARAMETERS")
end


### Check Auxiliary Parameters
params_ok = check_auxparams(hshiftmax, vshiftmax, rmedmax,
        boxwid, nit, irelonorm, vzmodel_type)
if input_ok
    println("Auxiliary parameters are ok!")
else
    println("ERROR: FIX AUXILIARY PARAMETERS")
end

### Print Input Parameters

println("\nInput verified! Check results below:")
println("====================================")
@printf("Input files:\n")
println("> Event list: ", inpD["fin_evlist"])
println("> Event format: ", inpD["evlist_fmt"])
println("> Station list: ", inpD["fin_stlist"])
println("> Station format: ", inpD["stlist_fmt"])
println("> Xcor dataset: ", inpD["fin_xcordat"])
println("> Xcor format: ", inpD["xcordat_fmt"])
println("> Tdif format: ", inpD["tdif_fmt"])
println("> Velocity model: ", inpD["fin_vzmdl"])
print("Travel-time table depths (min, max, space): ")
@printf("%.2f %.2f %.2f\n",inpD["tt_dep0"],
    inpD["tt_dep1"],inpD["tt_ddep"])
print("Travel-time table ranges (min, max, space): ")
@printf("%.2f %.2f %.2f\n",inpD["tt_del0"],
    inpD["tt_del1"],inpD["tt_ddel"])
@printf("GrowClust parameters:\n")
@printf("> rmin, delmax, rmsmax: %.3f %.1f %.3f\n",
    inpD["rmin"],inpD["delmax"],inpD["rmsmax"])
@printf("> rpsavgmin, rmincut, ngoodmin, iponly: % .3f %.3f %d %d\n",
    inpD["rpsavgmin"],inpD["rmincut"],inpD["ngoodmin"],inpD["iponly"])
@printf("> nboot, nbranch_min: %d %d\n",
    inpD["nboot"],inpD["nbranch_min"])
@printf("Output files:\n")
println("> Relocated catalog: ", inpD["fout_cat"])
println("> Cluster file: ", inpD["fout_clust"])
println("> Bootstrap file: ", inpD["fout_boot"])
println("> Log file: ", inpD["fout_log"])

### Read Catalog
println("\nReading event list:")
@time qdf = read_evlist(inpD["fin_evlist"],inpD["evlist_fmt"])
qid2qnum = Dict(zip(qdf.qid,qdf.qix))# maps event id to serial number
show(qdf)
println()

### Read Stations
print("\nReading station list")
@time sdf = read_stlist(inpD["fin_stlist"],inpD["stlist_fmt"])
show(sdf)
println()

### Read Xcor Data
println("\nReading xcor data")
@time xdf = read_xcordata(inpD,qdf[!,[:qix,:qid,:qlat,:qlon]],sdf[!,[:sta,:slat,:slon]])
show(xdf)
println()

###

### Read in velocity model
println("\nReading velocity model..")
# read in
z_s0, alpha_s0, beta_s0 = read_vzmodel(
    inpD["fin_vzmdl"],vpvs=inpD["vpvs_factor"])
for ii in 1:length(z_s0) # print out
    @printf("%5.2fkm: %6.4f %6.4f\n",z_s0[ii],alpha_s0[ii],beta_s0[ii])
end

### Find Moho depth in model, print results
println("\nMoho depths:")
imoho = find_moho(z_s0,alpha_s0)
@printf("%.2f %.3f %.3f\n",
    z_s0[imoho], alpha_s0[imoho], beta_s0[imoho])
println("Moho slownesses:")
@printf("%.4f %.4f\n",1.0/alpha_s0[imoho], 1.0/beta_s0[imoho])

### Interpolate VZ Model to finer grid
println("Interpolation:") # interpolate
z_s, alpha_s, beta_s = interp_vzmodel(
    z_s0, alpha_s0, beta_s0; itp_dz=inpD["itp_dz"])
for ii in 1:length(z_s) # print out
    @printf("%5.2fkm: %6.4f %6.4f\n",z_s[ii],alpha_s[ii],beta_s[ii])
end

### Run Earth-flattening Codes
erad = max(6371., z_s[end] + 0.1) # flatten
z, alpha = eflatten(z_s, alpha_s, erad=erad)
z, beta = eflatten(z_s, beta_s, erad=erad)


### Define slowness arrays
npts = length(z_s)
slow = zeros(Float64,npts,2)
slow[:,1] .= 1.0./alpha
slow[:,2] .= ifelse.(beta.>0.0,1.0./beta,1.0./alpha)

### Main Ray Tracing Loop
println("\nRay tracing to assemble travel time tables:")

# define grids
qdeptab = collect(range(inpD["tt_dep0"],inpD["tt_dep1"],step=inpD["tt_ddep"]))
sdeltab = collect(range(inpD["tt_del0"],inpD["tt_del1"],step=inpD["tt_ddel"]))

# loop over phases
phases = [1,2]
plongcuts = [inpD["plongcutP"],inpD["plongcutS"]]
ttoutfiles = [inpD["fout_pTT"], inpD["fout_sTT"]]
total_time = @elapsed for iphase in phases

    # Print results
    print("Working on Phase #")
    println(iphase)
    
    # Ray tracing: compute offset and travel time to different depths
    println("Tracing rays...")
    ptab, qdepxcor, qdeptcor, qdepucor, del2W, tt2W = trace_rays(
        iphase,z_s,z,slow,qdeptab,inpD["itp_dz"])
    println("Done.")
    
    # Make table of first arrivals and take of angles
    println("Compiling travel time table of first arrivals...")
    TT, AA = first_arrivals(vzmodel_type, plongcuts[iphase], qdeptab, sdeltab, 
                        slow[1,iphase], ptab, qdepxcor, qdeptcor, qdepucor, del2W, tt2W)
    println("Done.")
    println(sum(isnan.(TT)))
    
    # Write output files
    println("Writing output files...")
    write_table(ttoutfiles[iphase],inpD["fin_vzmdl"],iphase,vzmodel_type,
                        TT,qdeptab, sdeltab,ptab)
    println("Done.")
end
@printf("\nElapsed seconds: %.2f",total_time)
println()

### Test Interpolant objects

# instantiate interpolants
println("Testing smtrace interpolation (P and S waves):")
const pTT = smtrace_table(inpD["fout_pTT"],shallowmode,Float64)
const sTT = smtrace_table(inpD["fout_sTT"],shallowmode,Float64)
const ttTABs = [pTT,sTT]

# define test distances, depth
if vzmodel_type == 1
    test_dists = collect(range(0.0,80.0,step=5.0))
    npts = length(test_dists)
    test_depth = 10.0
else
    test_dists = collect(range(0.0,100.0,step=5.0))
    npts = length(test_dists)
    test_depth = 33.0
end
    
# calculate travel times
test_ttP = pTT(test_dists, test_depth)
test_ttS = sTT(test_dists, test_depth)

# print results
println(typeof(test_ttP[1]), " ", typeof(test_ttS[1]))
for ii in 1:npts
    if vzmodel_type == 1
        @printf("Distance: %4.0fkm, Depth: %3.0fkm",
            test_dists[ii],test_depth)
        @printf(" --> P and S travel times: %5.2fs, %5.2fs\n",
            test_ttP[ii],test_ttS[ii])
    else
        @printf("Distance: %4.0fdeg, Depth: %3.0fkm",
            test_dists[ii],test_depth)
        @printf(" --> P and S travel times: %5.3fmin, %5.3fmin\n",
            test_ttP[ii],test_ttS[ii])
    end
end
#exit()


#### Validate Event Depths and Travel Time Tables; Datum Setup

println("\nChecking event depths")

# print event depths
const min_qdep = minimum(qdf.qdep)
const max_qdep = maximum(qdf.qdep)
@printf("min and max event depth: %.3fkm %.3fkm\n",min_qdep,max_qdep)

# print table depths
@printf("min and max table depth: %.3fkm %.3fkm\n",inpD["tt_dep0"],inpD["tt_dep1"])

# implement warrnings and checks
if (min_qdep < inpD["tt_dep0"])
    println("WARNING: min event depth < min table depth")
end
if (min_qdep < z_s0[1])
    println("WARNING: min event depth < min vzmodel depth")
    #exit() # allow this, but warn user (should be ok if depth is near 0)
end
if (max_qdep > inpD["tt_dep1"]) # note tt_dep1 is >= vzmax
    println("ERROR: max event depth > max table / velocity model depth")
    exit()
end
if (inpD["tt_dep0"] < z_s0[1]) # for robustness, check this as well
    println("ERROR: min table depth < min vzmodel depth")
    exit()
end

# station elevations
if inpD["stlist_fmt"] == 2
    min_selev = minimum(sdf.selev)
    max_selev = maximum(sdf.selev)
    mean_selev = mean(sdf.selev)
    @printf("station elevation (min,mean,max): %.1fm %.1fm %.1fm\n",
        min_selev,mean_selev,max_selev)
    const datum = mean_selev/1000.0
else
    const datum = 0.0
end
@printf("assumed surface datum: %.1fkm\n",datum)

# check xcor data
println("\nValidating xcor data...")
nbad = sum(abs.(xdf.tdif).>tdifmax)
if nbad > 0
    println("Error: bad input differential times, some larger than tdifmax=$tdifmax")
    println("Fix input file or adjust tdifmax parameter.")
    ibad = abs.(xdf.tdif).>tdifmax
    show(xdf[ibad,:])
    exit()
end
nbad = sum(xdf.sdist.>inpD["tt_del1"])
if nbad > 0
    println("Error: bad input xcor data, stations further than travel time table allows")
    println("Fix input file or adjust travel time table parameter.")
    ibad = xdf.sdist.>inpD["tt_del1"]
    show(xdf[ibad,:])
    exit()
end
println("Done.")

############# Main Clustering Loop: Including Bootstrapping ##############

# loading packages @everywhere
@everywhere using Printf
@everywhere using DataFrames
@everywhere using StatsKit
@everywhere using Random
@everywhere using Dates
@everywhere using Interpolations
@everywhere using GrowClust1D
@everywhere using SharedArrays

# shared parameters - is there a better way?
@everywhere nboot = $(inpD["nboot"])
@everywhere nit = $nit
@everywhere boxwid = $boxwid
@everywhere degkm = $degkm
@everywhere irelonorm = $irelonorm
@everywhere rmsmax = $(inpD["rmsmax"])
@everywhere ngoodmin = $(inpD["ngoodmin"])
@everywhere rmedmax = $rmedmax
@everywhere distmax = $distmax
@everywhere distmax2 = $distmax2
@everywhere hshiftmax = $hshiftmax
@everywhere vshiftmax= $vshiftmax
@everywhere torgdifmax = $torgdifmax
@everywhere nupdate = $nupdate
@everywhere maxlink = $maxlink
@everywhere ttTABs = $ttTABs
@everywhere nq = $(nrow(qdf))
@everywhere iseed = $iseed

# initial locations
@everywhere qlat = $(qdf.qlat)
@everywhere qlon = $(qdf.qlon)
@everywhere qdep = $(qdf.qdep .+ datum)

# Setup bootstrapping matrices
blatM = SharedArray(repeat(qlat,1,nboot+1))
blonM = SharedArray(repeat(qlon,1,nboot+1))
bdepM = SharedArray(repeat(qdep,1,nboot+1))
borgM = SharedArray(zeros(Float64,(nq,nboot+1)))
bnbM = SharedArray(zeros(Int64,(nq,nboot+1)))
bcidM = SharedArray(repeat(Vector{Int32}(1:nq),1,nboot+1))
bnpairM = SharedArray(zeros(Int64,nboot+1))

# base xcor dataframe to sample from
xdf00 = select(xdf,[:qix1,:qix2,:slat,:slon,:tdif,:iphase,:igood])
xdf00[!,:gxcor] = ifelse.(xdf.igood.>0,xdf.rxcor,Float32(0.0)) # xcor with bad values zeroed
@everywhere xdf00 = $(xdf00)

# sampling vector
if nrow(xdf00)<typemax(Int32)
    @everywhere nxc = $(Int32(nrow(xdf00)))
    @everywhere ixc = Vector{Int32}(1:nxc)
else
    @everywhere nxc = $(nrow(xdf00))
    @everywhere ixc = Vector{Int64}(1:nxc)
end


#### loop over each bootstrapping iteration
println("\n\n\nStarting relocation estimates, workers=",workers())
@time @sync @distributed for ib in 0:nboot # need to call @sync to ensure all workers finish

    # log thread id
    @printf("Starting bootstrap iteration: %d/%d\n",ib,nboot)
    
    # timer for this thread
    wc = @elapsed begin
    
    # bootstrapping: resample data before run
    println("Initializing xcorr data and event pairs...")
    
    # bootstrapping: resample data before run
    Random.seed!(iseed + ib) # different for each run
    wc2 = @elapsed begin
    if ib > 0 # sample with replacement from original xcorr array
        isamp = sort(sample(ixc,nxc,replace=true)) # sorted to keep evpairs together
        rxdf = xdf00[isamp,:]
    else     # do not resample on first run
        rxdf = xdf00
    end

    # assign xcorr index array
    if nxc < typemax(Int32)
        rxdf[!,:ixx] = Vector{Int32}(1:nxc)
    else
        rxdf[!,:ixx] = Vector{Int64}(1:nxc)
    end

    # compile event pair arrays
    bpdf = combine(groupby(rxdf[!,Not([:slat,:slon,:iphase,:tdif])],[:qix1,:qix2]),
        :gxcor=>sum=>:rfactor,:ixx=>first=>:ix1,:ixx=>last=>:ix2,:igood=>sum=>:ngood)
    
    # keep only good pairs, and sort
    bpdf = bpdf[bpdf.ngood.>=ngoodmin,[:qix1,:qix2,:rfactor,:ix1,:ix2]]
    sort!(bpdf,:rfactor,rev=true)
    end # ends elapsed time for setup
    print("Done with initialization, elapsed time = $wc2\n")
    
    # run clustering
    brlats, brlons, brdeps, brorgs, brcids, bnbs = clustertree(
        bpdf.qix1, bpdf.qix2, bpdf.ix1, bpdf.ix2, 
        rxdf.tdif, rxdf.slat, rxdf.slon, rxdf.iphase,
        qlat, qlon, qdep, ttTABs, nit, boxwid, degkm, irelonorm,
        rmsmax,rmedmax,distmax,distmax2,
        hshiftmax,vshiftmax,torgdifmax,nupdate,maxlink)
    
    # save output to Shared Array
    if ib > 0
        blatM[:,ib] .= brlats
        blonM[:,ib] .= brlons
        bdepM[:,ib] .= brdeps
        borgM[:,ib] .= brorgs
        bnbM[:,ib] .= bnbs
        bcidM[:,ib] .= brcids
        bnpairM[ib] = nrow(bpdf)
    else
        blatM[:,nboot+1] .= brlats
        blonM[:,nboot+1] .= brlons
        bdepM[:,nboot+1] .= brdeps
        borgM[:,nboot+1] .= brorgs
        bnbM[:,nboot+1] .= bnbs
        bcidM[:,nboot+1] .= brcids
        bnpairM[nboot+1] = nrow(bpdf)
    end
        
    # completion
    end # ends the wall clock
    @printf("Completed bootstrap iteration: %d/%d, wall clock = %.1fs.\n",ib,nboot,wc)

end

################################################################

### Extract event-based output arrays
nq = Int32(nrow(qdf))
revids = qdf[:,:qid]
rlats = blatM[:,nboot+1]
rlons = blonM[:,nboot+1]
rdeps = bdepM[:,nboot+1]
rorgs = borgM[:,nboot+1]
rcids = bcidM[:,nboot+1]
npair = bnpairM[nboot+1]

### Finalize Clustering Trees ###

println("\nFinalizing clustering trees...")

# temporary dataframe with event and cluster number
tdf = DataFrame("enum"=>qdf.qix,"cnum"=>rcids)

# compute nbranch
transform!(groupby(tdf, :cnum), nrow => :nb)

# assign cluster ids, largest clusters 1st
sort!(tdf,:nb,rev=true)
tdf[!,:cid] .= 0
cc = 0
for sdf in groupby(tdf,:cnum)
    global cc+=1
    tdf[parentindices(sdf)[1],:cid] .= cc
end

# put this back in right order
sort!(tdf,:enum)
nclust = maximum(tdf.cid) 

# update cluster ids for all events
rcids00 = copy(rcids)
rnbranch = tdf.nb
rcids = tdf.cid

# finalize t-arrays
tlons, tlats = zeros(nclust), zeros(nclust) 
tdeps, torgs = zeros(nclust), zeros(nclust)
tnbranch = zeros(Int64,nclust)
for icc in 1:nclust
    idx = findall(rcids.==icc)
    tlons[icc] = mean(rlons[idx])
    tlats[icc] = mean(rlats[idx])
    tdeps[icc] = mean(rdeps[idx]) # note, datum-shifted
    torgs[icc] = mean(rorgs[idx])
    tnbranch[icc] = length(idx)
end

# compute cluster counts
ntree2 = sum(tnbranch.>=2)
ntree5 = sum(tnbranch.>=5)
ntree10 = sum(tnbranch.>=10)
ntree20 = sum(tnbranch.>=20)
ntree50 = sum(tnbranch.>=50)
ntree100 = sum(tnbranch.>=100)

### Relocated DataFrame
# (no datum adjustment yet, done on output)

# initialize
println("\nFinalizing relocated dataset:")
rdf = DataFrame("enum"=>qdf.qix,"evid"=>revids,
    "rlat"=>rlats,"rlon"=>rlons,"rdep"=>rdeps,"rtim"=>rorgs,
    "rcid"=>rcids,"rnb"=>rnbranch)

# compute relocated
nreloc = sum(rdf.rnb .> 1)
println("Relocated: ",nreloc)

# adjusted otime
rdf[!,:rot] = qdf[!,:qotime] .+ Nanosecond.(round.(rorgs*1.0e9))
show(rdf)
println()


### Compute Misfits - w/ otime adjustment ###

println("\nComputing misfits...")

# select columns
resdf = select(xdf,[:qid1,:qid2,:tdif,:iphase,:slat,:slon,:igood])
#resdf = resdf[resdf[!,:igood].==1,:] # f90 codes do this, not sure why

# merge with event data, renaming columns to specify 1/2
resdf = innerjoin(resdf,rdf[!,[:enum,:evid,:rlat,:rlon,:rdep,:rtim]],on=:qid1=>:evid)
DataFrames.rename!(resdf,:enum=>:qnum1,:rlat=>:qlat1,:rlon=>:qlon1,:rdep=>:qdep1,:rtim=>:qtim1)
resdf = innerjoin(resdf,rdf[!,[:enum,:evid,:rlat,:rlon,:rdep,:rtim]],on=:qid2=>:evid)
DataFrames.rename!(resdf,:enum=>:qnum2,:rlat=>:qlat2,:rlon=>:qlon2,:rdep=>:qdep2,:rtim=>:qtim2)

# compute source station distance
sdist1 = map_distance(resdf[!,:qlat1],resdf[!,:qlon1],resdf[!,:slat],resdf[!,:slon])
sdist2 = map_distance(resdf[!,:qlat2],resdf[!,:qlon2],resdf[!,:slat],resdf[!,:slon])

# compute predicted travel times
resdf[!,:pdif] = ifelse.(resdf[!,:iphase].==1,
    pTT.(sdist2,resdf[!,:qdep2]).-pTT.(sdist1,resdf[!,:qdep1]),
    sTT.(sdist2,resdf[!,:qdep2]).-sTT.(sdist1,resdf[!,:qdep1])) .+ 
    (resdf[!,:qtim2].-resdf[!,:qtim1]) # otime adjustment (add here or subtract from tdif)

# P vs S
ipp = (resdf[!,:iphase].==1)
npp = sum(ipp)
iss = .!ipp
nss = sum(iss)

# RMS
rmsP = evalrms(resdf[ipp,:tdif].-resdf[ipp,:pdif])
rmsS = evalrms(resdf[iss,:tdif].-resdf[iss,:pdif])

# Mean signed residuals
msresP = mean(resdf[ipp,:tdif].-resdf[ipp,:pdif])
msresS = mean(resdf[iss,:tdif].-resdf[iss,:pdif])

# print
println("P-wave RMS: $rmsP")
println("Mean signed residual: $msresP")
println("Phases used: $npp")
println("S-wave RMS: $rmsS")
println("Mean signed residual: $msresS")
println("Phases used: $nss")

# show
select!(resdf,[:qnum1,:qnum2,:qid1,:qid2,:tdif,:iphase,:pdif])

### Compute Event-based Stats

# arrays to store stats for each event
qnpair = zeros(Int64,nq)
qndiffP, qndiffS = zeros(Int64,nq), zeros(Int64,nq)
qsseP, qsseS = zeros(nq), zeros(nq)

# event 1 in pair, P-wave
for subdf in groupby(resdf[ipp,:],:qnum1)
    qix=subdf[1,:qnum1]
    qndiffP[qix] += nrow(subdf)
    qsseP[qix] += sum((subdf[!,:tdif].-subdf[!,:pdif]).^2)
end

# event 1 in pair, S-wave
for subdf in groupby(resdf[iss,:],:qnum1)
    qix=subdf[1,:qnum1]
    qndiffS[qix] += nrow(subdf)
    qsseS[qix] += sum((subdf[!,:tdif].-subdf[!,:pdif]).^2)
end

# event 2 in pair, P-wave
for subdf in groupby(resdf[ipp,:],:qnum2)
    qix=subdf[1,:qnum2]
    qndiffP[qix] += nrow(subdf)
    qsseP[qix] += sum((subdf[!,:tdif].-subdf[!,:pdif]).^2)
end

# event 2 in pair, S-wave
for subdf in groupby(resdf[iss,:],:qnum2)
    qix=subdf[1,:qnum2]
    qndiffS[qix] += nrow(subdf)
    qsseS[qix] += sum((subdf[!,:tdif].-subdf[!,:pdif]).^2)
end

# event pairs
for subdf in groupby(resdf,[:qnum1,:qnum2])
    qix1, qix2 = values(subdf[1,[:qnum1,:qnum2]])
    qnpair[qix1] += 1
    qnpair[qix2] += 1
end

# combined data --> RMS values
qrmsP = ifelse.(qndiffP.>0,sqrt.(qsseP./qndiffP),NaN64)
qrmsS = ifelse.(qndiffS.>0,sqrt.(qsseS./qndiffS),NaN64)

### Compute bootstrap statistics ###

# pre-allocate: defaults are NaN for errors, 0 for nb arrays
boot_stdH, boot_madH = fill(NaN64,nq), fill(NaN64,nq)
boot_stdZ, boot_madZ = fill(NaN64,nq), fill(NaN64,nq)
boot_stdT, boot_madT = fill(NaN64,nq), fill(NaN64,nq)
boot_nbL, boot_nbM, boot_nbH = zeros(Int64,nq), zeros(Float64,nq), zeros(Int64,nq)

# loop over events
if inpD["nboot"] > 1 # need to have run bootstrapping


    println("\nComputing bootstrap statistics...")
    for ii = 1:nq

        # nbranch statistics
        boot_nbL[ii] = minimum(bnbM[ii,:])
        boot_nbM[ii] = mean(bnbM[ii,:])
        boot_nbH[ii] = maximum(bnbM[ii,:])
        
        # only if relocated
        if rdf[ii,:rnb]>1

            # standard errors
            xstd = degkm*cosd(rdf[ii,:rlat])*std(blonM[ii,:])
            ystd = degkm*std(blatM[ii,:])
            boot_stdH[ii]=sqrt(xstd^2+ystd^2)
            boot_stdZ[ii]=std(bdepM[ii,:])
            boot_stdT[ii]=std(borgM[ii,:])

            # MAD statistics (no division by 1/quantile(Normal(), 3/4) ≈ 1.4826)
            xmad = degkm*cosd(rdf[ii,:rlat])*mad(blonM[ii,:],normalize=false)
            ymad = degkm*mad(blatM[ii,:],normalize=false)
            boot_madH[ii]=sqrt(xmad^2+ymad^2)
            boot_madZ[ii]=mad(bdepM[ii,:],normalize=false)
            boot_madT[ii]=mad(borgM[ii,:],normalize=false)

        end
    end
end
        
#########################################################

### Write Output File: Catalog

println("\nWriting output catalog: ", inpD["fout_cat"])

# open output file
fcat = open(inpD["fout_cat"],"w")

# loop over events
for ii = 1:nq
    
    # print out origin time, relocated position, magnitude
    dateS = Dates.format(rdf[ii,:rot],"YYYY mm dd HH MM SS.sss")
    @printf(fcat,"%s %9d %9.5f %10.5f %7.3f %5.2f ",
        dateS,rdf[ii,:evid],rdf[ii,:rlat],rdf[ii,:rlon],
        rdf[ii,:rdep]-datum,qdf[ii,:qmag]) # note, shifting datum back
    
    # print out cluster number and fits
    @printf(fcat,"%7d %7d %7d %5d %5d %5d %5.2f %5.2f ",
        rdf[ii,:enum],rdf[ii,:rcid],rdf[ii,:rnb],qnpair[ii],
        qndiffP[ii],qndiffS[ii],qrmsP[ii],qrmsS[ii])
    
    # print out uncertanties and catalog locations
    @printf(fcat,"%7.3f %7.3f %7.3f %9.5f %10.5f %7.3f\n",
        boot_madH[ii],boot_madZ[ii],boot_madT[ii],
        qdf[ii,:qlat],qdf[ii,:qlon],qdf[ii,:qdep])
end

# close file
close(fcat)



### Write Output File: Cluster

if !(inpD["fout_clust"] in ["none","NONE"])

    println("\nWriting output clusters")

    # open output file
    fcc = open(inpD["fout_clust"],"w")

    # loop over all clusters to output (only makes sense for n>=2)
    for cc in findall(tnbranch .>= min(2,inpD["nbranch_min"]))
        
        # write cluster info
        @printf(fcc,"%8d %7d %9.5f %10.5f %7.3f %7.3f\n",
            cc,tnbranch[cc],tlats[cc],tlons[cc],tdeps[cc]-datum,torgs[cc]) # datum-shift
        
        # write info for all events
        for ii in findall(rcids.==cc)
        
            # print out clustering, event info, mag, otime
            dateS = Dates.format(rdf[ii,:rot],"YYYY mm dd HH MM SS.sss")
            @printf(fcc,"%8d %8d %9d %5.2f %s ",cc,ii,rdf[ii,:evid],qdf[ii,:qmag],dateS)
            
            # print event location and position w/in cluster
            qdx = degkm*(rdf[ii,:rlon]-tlons[cc])*cosd(tlats[cc])
            qdy = degkm*(rdf[ii,:rlat]-tlats[cc])
            qdz = rdf[ii,:rdep]-tdeps[cc] # neither are datum-adjusted in array
            @printf(fcc,"%9.5f %10.5f %7.3f %9.4f %9.4f %9.4f ",
                rdf[ii,:rlat],rdf[ii,:rlon],rdf[ii,:rdep]-datum,qdx,qdy,qdz)

            # print out uncertanties  and catalog locations
            @printf(fcc,"%7.3f %7.3f %7.3f %9.5f %10.5f %7.3f\n",
                boot_madH[ii],boot_madZ[ii],boot_madT[ii],
                qdf[ii,:qlat],qdf[ii,:qlon],qdf[ii,:qdep])
        end
        
    end


    # close file 
    close(fcc)

end

### Write Output File: Bootstrapping

# only if requested
if (inpD["nboot"] > 1)&(!(inpD["fout_boot"] in ["none","NONE"]))
    
    println("\nWriting output bootstrapping")
    
    # open file
    fbb = open(inpD["fout_boot"],"w")
    
    # file header
    @printf(fbb,"%8d %5d\n",nq,inpD["nboot"])
    
    
    # loop over events
    for ii = 1:nq
        
        # event header
        @printf(fbb,"%8d %9d %9.5f %10.5f %7.3f ",
            rdf[ii,:enum],rdf[ii,:evid],rdf[ii,:rlat],rdf[ii,:rlon],rdf[ii,:rdep]-datum)
        @printf(fbb,"%7d %7d %6.2f %7d %7d %9.5f %10.5f %7.3f ",
            rdf[ii,:rcid],rdf[ii,:rnb],boot_nbM[ii],boot_nbL[ii],boot_nbH[ii],
            qdf[ii,:qlat],qdf[ii,:qlon],qdf[ii,:qdep])
        @printf(fbb,"%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
            boot_madH[ii],boot_madZ[ii],boot_madT[ii],
            boot_stdH[ii],boot_stdZ[ii],boot_stdT[ii])
        
        # report each bootstrap (f90 version reports CID but this is dumb?)
        for ib = 1:inpD["nboot"]
            @printf(fbb,"%9.5f %10.5f %7.3f %7d\n",
                blatM[ii,ib],blonM[ii,ib],bdepM[ii,ib],bnbM[ii,ib])
        end
        
    end
    
    # close file
    close(fbb) 
    
end


### Write Output File: Log / Statistics

println("\nWriting run log")

# open log file
flog = open(inpD["fout_log"],"w")

# Run Parameters
@printf(flog,  "************************ Input files ************************\n")
@printf(flog,  "     control file:   %s\n", infile_ctl)
@printf(flog,  "       event list:   %s\n", inpD["fin_evlist"])
@printf(flog,  "     station list:   %s\n", inpD["fin_stlist"])
@printf(flog,  "     xcordat file:   %s\n", inpD["fin_xcordat"])
@printf(flog,  "     velocity mdl:   %s\n", inpD["fin_vzmdl"])
@printf(flog, "\n")
@printf(flog,  "********************** Travel Time Tables *********************\n")
@printf(flog,  "          P-phase:   %s\n", inpD["fout_pTT"])
@printf(flog,  "          S-phase:   %s\n", inpD["fout_sTT"])
@printf(flog, "\n")
@printf(flog,  "************************ Output files *************************\n")
@printf(flog,  "     catalog file:   %s\n", inpD["fout_cat"])
@printf(flog,  "     cluster file:   %s\n", inpD["fout_clust"])
@printf(flog,  "         log file:   %s\n", inpD["fout_log"])
@printf(flog,  "   bootstrap file:   %s\n", inpD["fout_boot"])
@printf(flog, "\n")
@printf(flog,  "****************** GROWCLUST Run Parameters *******************\n")
@printf(flog, "%56s %6.2f\n", " (min rxcor value for evpair similarity coeff.): rmin =", inpD["rmin"])
@printf(flog, "%56s %6.1f\n", " (max sta. dist for evpair similarity coeff.): delmax =", inpD["delmax"])
@printf(flog, "%56s %6.2f\n", " (max rms residual to join clusters): rmsmax =", inpD["rmsmax"])
@printf(flog, "%56s %6d\n" , " (num. bootstrap uncertainty iterations): nboot =", inpD["nboot"])
@printf(flog, "\n")
@printf(flog,  "****************** Auxiliary Run Parameters *******************\n")
@printf(flog, "%56s %6.2f\n", " max catalog dist to join clusters: ", distmax)
@printf(flog, "%56s %6.2f\n", " max relocated dist to join clusters: ", distmax2)
@printf(flog, "%56s %6.2f\n", " max permitted horizontal cluster shifts: ", hshiftmax)
@printf(flog, "%56s %6.2f\n", " max permitted vertical cluster shifts: ", vshiftmax)
@printf(flog, "%56s %6.2f\n", " max median absolute residual to join clusters: ", rmedmax)  
@printf(flog,  "***************************************************************\n")


# Run Summary with statistics
@printf(flog, "\n")
@printf(flog, "==================================================================\n")
@printf(flog, "==================================================================\n")
@printf(flog, "\n")
@printf(flog, "********************  GROWCLUST Run Summary  *********************\n")
#@printf(flog, "\n")
@printf(flog, "%55s %10d\n", "Number of catalog events: ", nq)
@printf(flog, "%55s %10d\n", "Number of relocated events: ", nreloc)
@printf(flog, "%55s %10d\n", "Number of input event pairs: ", npair)
@printf(flog, "%55s %10d\n", "Number of event pairs used: ", sum(qnpair)/2)
@printf(flog, "%55s %10d\n", "Number of xcor data used (total, P+S): ", npp + nss)
@printf(flog, "%55s %10d\n", "Number of xcor data used (P-phase): ", npp)
@printf(flog, "%55s %10.4f\n", "RMS differential time residual (P-phase): ", rmsP)
@printf(flog, "%55s %10.4f\n", "Mean (signed) differential time residual (P-phase): ", msresP)
@printf(flog, "%55s %10d\n", "Number of xcor data used (S-phase): ", nss)
@printf(flog, "%55s %10.4f\n", "RMS differential time residual (S-phase): ", rmsS)
@printf(flog, "%55s %10.4f\n", "Mean (signed) differential time residual (S-phase): ", msresP)
@printf(flog,  "\n")
@printf(flog, "%55s %9d\n", "Number of clusters with >=   2 events: ", ntree2)
@printf(flog, "%55s %9d\n", "Number of clusters with >=   5 events: ", ntree5)
@printf(flog, "%55s %9d\n", "Number of clusters with >=  10 events: ", ntree10)
@printf(flog, "%55s %9d\n", "Number of clusters with >=  20 events: ", ntree20)
@printf(flog, "%55s %9d\n", "Number of clusters with >=  50 events: ", ntree50)
@printf(flog, "%55s %9d\n", "Number of clusters with >= 100 events: ", ntree100)
   
# close this
close(flog)
