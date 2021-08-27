### Script for Running GrowClust: Julia

### Import packages

# Packages
using Printf
using DataFrames
using CSVFiles
using StatsKit
using Random
using Dates
using Interpolations

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
const nupdate = 100         # update progress every nupdate pairs - NEW
   
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

# define event-based output arrays
const nq = Int32(nrow(qdf))
revids = qdf[:,:qid]
rlats = qdf[:,:qlat]
rlons = qdf[:,:qlon]
rdeps = qdf[:,:qdep] .+ datum # datum-adjust
rorgs = zeros(Float64,nq) # origin time adjust
rcids = Vector{Int32}(1:nq) # initialize each event into one cluster

# Setup bootstrapping matrices
if inpD["nboot"] > 0
    blatM = repeat(qdf.qlat,1,inpD["nboot"])
    blonM = repeat(qdf.qlon,1,inpD["nboot"])
    bdepM = repeat(qdf.qdep,1,inpD["nboot"]) .+ datum
    borgM = zeros(Float64,(nq,inpD["nboot"]))
    bnbM = repeat(Vector{Int32}(1:nq),1,inpD["nboot"])
end

# base xcor dataframe to sample from
xdf00 = select(xdf,[:qix1,:qix2,:slat,:slon,:tdif,:iphase,:igood])
xdf00[!,:gxcor] = ifelse.(xdf.igood.>0,xdf.rxcor,Float32(0.0)) # xcor with bad values zeroed
#show(xdf00)

# sampling vector
if nrow(xdf00)<typemax(Int32)
    const nxc = Int32(nrow(xdf00))
    const ixc = Vector{Int32}(1:nxc)
else
    const nxc = nrow(xdf00)
    const ixc = Vector{Int64}(1:nxc)
end

# parameter constants
const rmsmax = inpD["rmsmax"]
const cdepmin = min(0.0,min_qdep)
const hshiftmaxD2 = (hshiftmax/degkm)^2 # in squared degrees
const distmax22 = distmax^2 # in squared km


# loop over each bootstrapping iteration
Random.seed!(iseed)
println("\n\n\nStarting relocation estimates, nthread=",Threads.nthreads())
@time Threads.@threads for ib in 0:inpD["nboot"]   
    
    # log thread id
    @printf("Thread %d: starting bootstrap iteration: %d/%d\n",
            Threads.threadid(),ib,inpD["nboot"])
    
    # timer for this thread
    wc = @elapsed begin
    
    # bootstrapping: resample data before run
    if Threads.threadid()==1
        println("Thread 1: Initializing xcorr data and event pairs...")
    end
    if ib > 0
        
        # sample from original xcorr array
        isamp = sort(sample(ixc,nxc,replace=true)) # sorted to keep evpairs together
        rxdf = xdf00[isamp,:]
        if nxc < typemax(Int32)
            rxdf[!,:ixx] = Vector{Int32}(1:nxc)
        else
            rxdf[!,:ixx] = Vector{Int64}(1:nxc)
        end

        # define event-based arrays
        brlats = blatM[:,ib]
        brlons = blonM[:,ib]
        brdeps = bdepM[:,ib]
        brorgs = borgM[:,ib]
        brevids = qdf[!,:qid]
        brcids = Vector{Int32}(1:nq) # initialize each event into one cluster

    # do not resample on first run
    else
        
        # use original xcorr array
        rxdf = xdf00
        if nxc < typemax(Int32)
            rxdf[!,:ixx] = Vector{Int32}(1:nxc)
        else
            rxdf[!,:ixx] = Vector{Int64}(1:nxc)
        end

        # precomputed event-based arrays
        brlats = rlats
        brlons = rlons
        brdeps = rdeps
        brorgs = rorgs
        brevids = revids
        brcids = rcids

    end

    # event pair arrays
    bpdf = combine(groupby(rxdf[!,Not([:slat,:slon,:iphase,:tdif])],[:qix1,:qix2]),
        :gxcor=>sum=>:rfactor,:ixx=>first=>:ix1,:ixx=>last=>:ix2,:igood=>sum=>:ngood)
    
    # keep only good pairs, and sort
    bpdf = bpdf[bpdf.ngood.>=inpD["ngoodmin"],[:qix1,:qix2,:rfactor,:ix1,:ix2]]
    sort!(bpdf,:rfactor,rev=true)
    bnpair = Int32(nrow(bpdf))
    #show(bpdf)
    
    
    # initialize clustering tree arrays
    btlats = copy(brlats)
    btlons = copy(brlons)
    btdeps = copy(brdeps)
    btorgs = copy(brorgs)
    btnbranch = ones(Int32,nq)

    # dictionary with cid => event indices (note cid=qix on initialization)
    cid2qixD = Dict(qix=>[qix] for qix in brcids)

    # dictionary with cid => event pair indices
    cid2pairD = Dict( qix => Vector{Int32}(findall(
        (bpdf.qix1.==qix).|(bpdf.qix2 .== qix))) for qix in brcids) # cid=qix on initialization

    # loop over event pairs
    @inbounds for ip in Vector{Int32}(1:bnpair)

        # Progress
        if ((mod(ip,nupdate)==0) & (Threads.threadid()==1))
            println("Thread 1: working on sorted pair: $ip/$bnpair")
        end
          
        # get event pair information
        qix1, qix2, rfactor = values(bpdf[ip,:])
        qc1, qc2 = brcids[qix1], brcids[qix2]
        if qc1 == qc2; continue; end # skip if in same cluster
        nb1, nb2 = btnbranch[qc1], btnbranch[qc2]

        # check to see if clusters are too far apart
        cdist22 = (btdeps[qc2]-btdeps[qc1])^2 + map_distance(
                btlats[qc1],btlons[qc1],btlats[qc2],btlons[qc2])^2
        if cdist22 > distmax22
            continue
        end

        # find event ids in each cluster
        c1qixs = cid2qixD[qc1] # possible speedup, not sure for this dataset  
        c2qixs = cid2qixD[qc2] # possible speedup, not sure for this dataset

        # only one link
        if (nb1==1)&(nb2==1)

            # only one link
            nlink = 1
            linx = [ip]

        # find all links between clusters
        else

            # # event pairs linking clusters
            # #   - pairs above are either processed and in same cluster
            # #     or caused issues when trying to link
            cpix1, cpix2 = cid2pairD[qc1],cid2pairD[qc2]
            @views linx = intersect(cpix1[cpix1.>=ip],cpix2[cpix2.>=ip])

            # keeping only best linx
            nlink = length(linx)
            if nlink > maxlink
                linx = partialsort(linx,1:maxlink)
            end 

        end

        # count number of picks
        npick = sum([bpdf[jj,:ix2]-bpdf[jj,:ix1]+1 for jj in linx])

        # extract locations relative to centroid
        dqlat1, dqlon1 = zeros(npick),zeros(npick)
        dqdep1, dqorg1 = zeros(npick),zeros(npick)
        dqlat2, dqlon2 = zeros(npick),zeros(npick)
        dqdep2, dqorg2 = zeros(npick),zeros(npick)
        phase12, slat12, slon12 = zeros(Int8,npick), zeros(npick), zeros(npick)
        tdif = zeros(npick)
        ix1 = 1 # start index for event pair in the arrays above
        @inbounds for ilink in linx
            pix1, pix2, jx1, jx2 = values(bpdf[ilink,[:qix1,:qix2,:ix1,:ix2]])
            ix2 = ix1 + (jx2 - jx1) # end-index in npick array
            if brcids[pix1]==qc1 # regular: event 1 in cluster 1, event 2 in cluster 2
                dqlat1[ix1:ix2] .= brlats[pix1]-btlats[qc1]
                dqlat2[ix1:ix2] .= brlats[pix2]-btlats[qc2]
                dqlon1[ix1:ix2] .= brlons[pix1]-btlons[qc1]
                dqlon2[ix1:ix2] .= brlons[pix2]-btlons[qc2]
                dqdep1[ix1:ix2] .= brdeps[pix1]-btdeps[qc1]
                dqdep2[ix1:ix2] .= brdeps[pix2]-btdeps[qc2]
                dqorg1[ix1:ix2] .= brorgs[pix1]-btorgs[qc1]
                dqorg2[ix1:ix2] .= brorgs[pix2]-btorgs[qc2]
                tdif[ix1:ix2] .= rxdf.tdif[jx1:jx2] # keep the tdif with same sign
            else # flipped: event 1 in cluster 2, event 2 in cluster 1
                dqlat1[ix1:ix2] .= brlats[pix2]-btlats[qc1]
                dqlat2[ix1:ix2] .= brlats[pix1]-btlats[qc2]
                dqlon1[ix1:ix2] .= brlons[pix2]-btlons[qc1]
                dqlon2[ix1:ix2] .= brlons[pix1]-btlons[qc2]
                dqdep1[ix1:ix2] .= brdeps[pix2]-btdeps[qc1]
                dqdep2[ix1:ix2] .= brdeps[pix1]-btdeps[qc2]
                dqorg1[ix1:ix2] .= brorgs[pix2]-btorgs[qc1]
                dqorg2[ix1:ix2] .= brorgs[pix1]-btorgs[qc2]
                tdif[ix1:ix2] .= -rxdf.tdif[jx1:jx2] # flip the tdif
            end
            slat12[ix1:ix2] .= rxdf.slat[jx1:jx2] # stays same
            slon12[ix1:ix2] .= rxdf.slon[jx1:jx2] # stays same
            phase12[ix1:ix2] .= rxdf.iphase[jx1:jx2] # stays same
            ix1 = ix2 + 1 # update start index in npick array
        end
        
        # unweighted cluster centroid
        clat0 = (btlats[qc1] + btlats[qc2]) / 2.0
        clon0 = (btlons[qc1] + btlons[qc2]) / 2.0
        cdep0 = (btdeps[qc1] + btdeps[qc2]) / 2.0

        # run difclust (relocation norms 1, 2, 3)
        if irelonorm == 1
            clat1, clon1, cdep1, clat2, clon2, cdep2, 
                cdist, torgdif, resid, rms, rmed, resol = difclust1(
                clat0,clon0,cdep0,tdif,phase12, slat12, slon12,
                dqlat1, dqlon1, dqdep1, dqorg1, dqlat2, dqlon2, dqdep2, dqorg2,
                ttTABs,boxwid,nit,degkm)
        elseif irelonorm == 2
            clat1, clon1, cdep1, clat2, clon2, cdep2, 
                cdist, torgdif, resid, rms, rmed, resol = difclust2(
                clat0,clon0,cdep0,tdif,phase12, slat12, slon12,
                dqlat1, dqlon1, dqdep1, dqorg1, dqlat2, dqlon2, dqdep2, dqorg2,
                ttTABs,boxwid,nit,degkm)
        else
            clat1, clon1, cdep1, clat2, clon2, cdep2, 
            cdist, torgdif, resid, rms, rmed, resol = difclust3(
                clat0,clon0,cdep0,tdif,phase12, slat12, slon12,
                dqlat1, dqlon1, dqdep1, dqorg1, dqlat2, dqlon2, dqdep2, dqorg2,
                ttTABs,boxwid,nit,degkm)
        end

        # careful with cluster mergers near surface
        if (min(cdep1,cdep2) < cdepmin)
            continue
        end

        # for robustness
        if abs(torgdif > torgdifmax)
            println("LARGE ORIGIN TIME CORRECTION: $torgdif")
            println("Likely xcor data or event list error!")
            qid1, qid2 = revids[qix1], revids[qix2]
            println("Current pair: $qid1 $qid2")
            show(bestdf)
            exit()
        end

        # reject cluster merger if rms or median absolute residual too large or relocated dist to far
        if ((rms > rmsmax) | (rmed > rmedmax) | (cdist > distmax2))
            continue
        end

        # fraction of events in each cluster
        fracC1 = Float64(nb1)/Float64(nb1+nb2)
        fracC2 = 1.0-fracC1

        # original centroid of combined cluster
        #creflat = cosd(clat0) # needed to scale dlon
        cxlat00 = btlats[qc1]*fracC1+btlats[qc2]*fracC2
        cxlon00 = btlons[qc1]*fracC1+btlons[qc2]*fracC2
        cxdep00 = btdeps[qc1]*fracC1+btdeps[qc2]*fracC2

        # new centroid of combined cluster
        cxlat11 = clat1*fracC1+clat2*fracC2
        cxlon11 = clon1*fracC1+clon2*fracC2
        cxdep11 = cdep1*fracC1+cdep2*fracC2

        # offset between two (b/c not all links used)
        dcxlat = cxlat11-cxlat00
        dcxlon = cxlon11-cxlon00
        dcxdep = cxdep11-cxdep00


        # check relative shift of cluster 1 (subtracting possible DC offset)
        # offsets for cluster 1: new loc - old loc
        qlat_off1 = clat1 - btlats[qc1]
        qlon_off1 = clon1 - btlons[qc1]
        qdep_off1 = cdep1 - btdeps[qc1]        
        if abs(qdep_off1-dcxdep) > vshiftmax
            continue
        elseif ((qlon_off1-dcxlon)*cosd(clat0))^2 + (qlat_off1-dcxlat)^2 > hshiftmaxD2 # in squared degrees
            continue
        end 

        # check relative shift of cluster 2 (subtracting possible DC offset)
        # offsets for cluster 2: new loc - old loc
        qlat_off2 = clat2 - btlats[qc2]
        qlon_off2 = clon2 - btlons[qc2]
        qdep_off2 = cdep2 - btdeps[qc2]
        if abs(qdep_off2 - dcxdep) > vshiftmax
            continue
        elseif ((qlon_off2-dcxlon)*cosd(clat0))^2 + (qlat_off2-dcxlat)^2 > hshiftmaxD2 # in squared degrees
            continue
        end 

        # ok, cluster merger is approved!
        #println("Approved!")

        # origin time updates
        @inbounds begin
        #qtim_off1 = btorgs[qc1] - torgdif/2.0 # symmetric shift "backward" if positive
        #qtim_off2 = btorgs[qc2] + torgdif/2.0 # symmetric shift "forward" if positive
        qtim_off1 = -torgdif/2.0 # symmetric shift "backward" if positive (torg always 0)
        qtim_off2 =  torgdif/2.0 # symmetric shift "forward" if positive (torg always 0)
        cxorg11 = fracC1*qtim_off1 + fracC2*qtim_off2

        # update locations: cluster 1
        brlats[c1qixs] .+= qlat_off1
        brlons[c1qixs] .+= qlon_off1
        brdeps[c1qixs] .+= qdep_off1
        brorgs[c1qixs] .+= qtim_off1

        # update locations: cluster 2
        brlats[c2qixs] .+= qlat_off2
        brlons[c2qixs] .+= qlon_off2
        brdeps[c2qixs] .+= qdep_off2
        brorgs[c2qixs] .+= qtim_off2

        # evacuate tree 2, assign events to tree 1
        brcids[c2qixs] .= qc1
        btnbranch[qc1] += nb2
        #btnbranch[qc2] = 0  # correct, but not needed

        end # end of @inbounds

        # merged cluster
        union!(cid2pairD[qc1],cid2pairD[qc2])
        union!(cid2qixD[qc1],c2qixs)
        

        # # align cluster and catalog centroids, set average time to zero
        # #  (need to do this b/c not all events used in difclust
        iclust = cid2qixD[qc1]
        @inbounds begin
        brorgs[iclust] .-= cxorg11 #mean(brorgs[iclust])
        brlats[iclust] .+= (cxlat00 - cxlat11)
        btlats[iclust] .= cxlat00
        brlons[iclust] .+= (cxlon00 - cxlon11)
        btlons[iclust] .= cxlon00
        brdeps[iclust] .+= (cxdep00 - cxdep11)
        btdeps[iclust] .= cxdep00
        #btorgs[iclust] .= 0.0 # always zero, never updated
        end
        
    end
    
    # save output
    if ib > 0
        blatM[:,ib] .= brlats
        blonM[:,ib] .= brlons
        bdepM[:,ib] .= brdeps
        borgM[:,ib] .= brorgs
        bnbM[:,ib] .= btnbranch[brcids]
    else
        rlats .= brlats
        rlons .= brlons
        rdeps .= brdeps
        rorgs .= brorgs
        rcids .= brcids
        global npair = bnpair
    end
        
    # completion
    end # ends the wall clock
    @printf("Thread %d: completed bootstrap iteration: %d/%d, wall clock = %.1fs.\n",
            Threads.threadid(),ib,inpD["nboot"],wc)
    println()
end

################################################################

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


### Compute Misfits - much better w/ otime adjustment ###

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


# loop over events
if inpD["nboot"] > 1 # need to have run bootstrapping


    println("\nComputing bootstrap statistics...")

    # pre-allocate: defaults are NaN for errors, 0 for nb arrays
    boot_stdH, boot_madH = fill(NaN64,nq), fill(NaN64,nq)
    boot_stdZ, boot_madZ = fill(NaN64,nq), fill(NaN64,nq)
    boot_stdT, boot_madT = fill(NaN64,nq), fill(NaN64,nq)
    boot_nbL, boot_nbM, boot_nbH = zeros(Int64,nq), zeros(Float64,nq), zeros(Int64,nq)


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