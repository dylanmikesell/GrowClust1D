### Simple functions to compute distances between points on a map ####
#  > Inputs: lat1, lon1, lat2, lon2 (float or vector)
#  < Returns: distance in km  (float or vector)

# scalar to scalar
function map_distance(lat1::Float64,lon1::Float64,
    lat2::Float64,lon2::Float64,degkm=111.1949266)
lat0 = 0.5*(lat2+lat1)
return degkm*sqrt( (lat2-lat1)^2 + (cosd(lat0)*(lon2-lon1))^2 )
end

# scalar to vector
function map_distance(lat1::Float64,lon1::Float64,
    lat2::Vector{Float64},lon2::Vector{Float64},degkm=111.1949266)
lat0 = 0.5*(mean(lat2)+lat1)
return degkm*sqrt.( (lat2.-lat1).^2 .+ (cosd(lat0)*(lon2.-lon1)).^2 )
end

# vector to scalar
function map_distance(lat1::Vector{Float64},lon1::Vector{Float64},
    lat2::Float64,lon2::Float64,degkm=111.1949266)
lat0 = 0.5*(mean(lat1)+lat2)
return degkm*sqrt.( (lat2.-lat1).^2 .+ (cosd(lat0)*(lon2.-lon1)).^2 )
end

# vector to vector
function map_distance(lat1::Vector{Float64},lon1::Vector{Float64},
    lat2::Vector{Float64},lon2::Vector{Float64},degkm=111.1949266)
lat0 = 0.5*mean(lat2.+lat1)
return degkm*sqrt.( (lat2.-lat1).^2 .+ (cosd(lat0)*(lon2.-lon1)).^2 )
end

### Robust Mean Based on Huber norm

function robomean(xx::Vector{Float64},xgap::Float64,nit::Int64)
    
    # get npts
    ndat = length(xx)
    
    # set initial weight
    xw = ones(Float64,ndat)
    
    # iterations loop
    xmean = 0.0
    for it = 1:nit
        
        # compute weighted mean
        xmean = sum(xw.*xx)/sum(xw)
        
        # compute weights for next iteration
        if it < nit # distance is xgap + |xx-xmean|, w = 1/d
            xw .= 1.0 ./ (xgap.+abs.(xx.-xmean))
        end
        
    end
    
    # compute final weighted misfit
    fit2 = sum(xw .* (xx.-xmean).^2)
    
    # return
    return xmean, fit2 
end

function robomean(xx::Vector{Float32},xgap::Float32,nit::Int64)
    
    # get npts
    ndat = length(xx)
    
    # set initial weight
    xw = ones(Float32,ndat)
    
    # iterations loop
    xmean = Float32(0.0)
    for it = 1:nit
        
        # compute weighted mean
        xmean = sum(xw.*xx)/sum(xw)
        
        # compute weights for next iteration
        if it < nit # distance is xgap + |xx-xmean|, w = 1/d
            xw .= Float32(1.0) ./ (xgap.+abs.(xx.-xmean))
        end
        
    end
    
    # compute final weighted misfit
    fit2 = sum(xw .* (xx.-xmean).^2)
    
    # return
    return xmean, fit2 
end

### Define rms function
function evalrms(data)
    return sqrt.(mean(data.^2))
end


# DIFCLUST performs relative relocation of two clusters of events
# (relative to the centroid of the cluster pair) using the the npr "best" 
# event pairs linking the clusters. (npr is not input explicitly,
# as the npick differential travel times, etc., include all observations 
# across these linking pairs). 
# This is analogous to DIFLOC, which does the same thing for two events, 
# relative to their centroid... As in DIFLOC, the method uses an iterative
# ("shrinking-box") grid search approach to obtain the best (L1) relative locations.
#----------
# Three versions: difclust1, difclust2, difclust3 --> for L1, L2, and Robomean norm.

# Inputs: qlat0  =  reference center point latitude
#         qlon0  =  reference center point longitude
#         qdep0  =  reference center point depth (km)
#         tdif   =  array (len=npick) of dif times, t2-t1 (s)
#         iph    =  array (len=npick) with phase index numbers (1 to 10) for tt data
#         slat   =  array (len=npick) with station latitudes
#         slon   =  array (len=npick) with station longitudes
#         qlat1  =  array (len=npick) of events in cluster1 latitude offsets from centroid
#         qlon1  =  array (len=npick) of events in cluster1 longitude offsets
#         qdep1  =  array (len=npick) of events in cluster1 depth offsets
#         qorg1  =  array (len=npick) of events in cluster1 time offsets
#         qlat2  =  array (len=npick) of events in cluster2 latitude offsets from centroid
#         qlon2  =  array (len=npick) of events in cluster2 longitude offsets
#         qdep2  =  array (len=npick) of events in cluster2 depth offsets
#         qorg2  =  array (len=npick) of events in cluster2 time offsets
#         ttTABs =  travel time tables for P and S phases
#         boxwid =  starting box width (km)
#         nit    =  number of iterations to perform
#         degkm  =  degrees to km conversion factor
# Returns: clat1  =  best-fitting latitude for first cluster centroid
#          clon1  =  best-fitting longitude for first cluster centroid
#          cdep1  =  best-fitting depth (km) of first cluster centroid
#          clat2  =  best-fitting latitude for second cluster centroid
#          clon2  =  best-fitting longitude for second cluster centroid
#          cdep2  =  best-fitting depth (km) of second cluster centroid
#          cdist  =  cluster separation distance (km)
#          torgdif=  origin time difference, i.e., t2-t1 median residual (>0 when 2 is later than 1)
#          resid  =  array (len=npick) of residuals (s) between observed tdif (tt) and predicted
#          rms    =  rms residual ( sqrt( sum(resid**2)/npick) )
#          rmed   =  median absolute value residual ( median( abs(resid) ) )
#          resol  =  nominal resolution (m) of final box

##### difclust1: L1 residual norm
function difclust1(qlat0::Float64,qlon0::Float64,qdep0::Float64,
    tdif::Vector{Float64},iph::Vector{Int8},
    slat::Vector{Float64},slon::Vector{Float64},
    qlat1::Vector{Float64},qlon1::Vector{Float64},
    qdep1::Vector{Float64},qorg1::Vector{Float64},
    qlat2::Vector{Float64},qlon2::Vector{Float64},
    qdep2::Vector{Float64},qorg2::Vector{Float64},
    ttTABs,boxwid::Float64,nit::Int64,degkm::Float64)

# initialize some variables
flatbest1, flonbest1, fdepbest1 = 0.0, 0.0, 0.0
flatbest2, flonbest2, fdepbest2 = 0.0, 0.0, 0.0
torgdif = 0.0

# extract npick
npick = length(tdif)
resid = zeros(npick)

# initialize box
dlat0 = 0.0
dlon0 = 0.0
ddep0 = 0.0
dlat = 0.5*boxwid/degkm
cosqlat = cosd(qlat0)
dlon = dlat/cosqlat

# careful near surface (added 07/2018)
zboxwid = boxwid
if (qdep0 < 0.0)
    zboxwid = min(1.0, boxwid)
elseif (qdep0 < boxwid/2.0)
    zboxwid = max(1.0,2.0*qdep0)
end
ddep = 0.5*zboxwid

# start iterations
for it=1:nit
    
    # grid search over box for best fit
    fitbest = 1.0e20
    tbest = 0.0
    
    # get trial locations: f1, f2 (3x3 grid for this box iteration)
    for iy = -1.0:1.0
        flat1 = qlat0 + dlat0 + dlat*iy
        flat2 = qlat0 - dlat0 - dlat*iy
        for ix = -1.0:1.0
            flon1 = qlon0 + dlon0 + dlon*ix
            flon2 = qlon0 - dlon0 - dlon*ix
            for iz = -1.0:1.0
                fdep1 = qdep0 + ddep0 + ddep*iz
                fdep2 = qdep0 - ddep0 - ddep*iz
                
                # compute predicted travel time and residuals w/observed
                sdist1 = map_distance(qlat1.+flat1,qlon1.+flon1,slat,slon)
                sdist2 = map_distance(qlat2.+flat2,qlon2.+flon2,slat,slon)                
                @inbounds for ii=1:npick
                    tt1 = ttTABs[iph[ii]](sdist1[ii],qdep1[ii]+fdep1)
                    tt2 = ttTABs[iph[ii]](sdist2[ii],qdep2[ii]+fdep2)
                    pdif = (tt2 + qorg2[ii]) - (tt1 + qorg1[ii])
                    resid[ii] = tdif[ii] - pdif # accounts for otime adjustment
                end
                
                # compute fit: L1 norm
                residval = median(resid)
                fit = sum(abs.(resid.-residval))
                
                # update best fit
                if fit < fitbest
                    fitbest = fit
                    flatbest1 = flat1
                    flonbest1 = flon1
                    fdepbest1 = fdep1
                    flatbest2 = flat2
                    flonbest2 = flon2
                    fdepbest2 = fdep2
                    tbest = residval
                end
                
            end
        end
    end # end of grid search
    
    # update best position
    dlat0 = flatbest1 - qlat0
    dlon0 = flonbest1 - qlon0
    ddep0 = fdepbest1 - qdep0
    torgdif = tbest # otime difference: residual between observed and predicted tdif

    # shrink box by 2/3 each iteration
    dlat *= 2.0/3.0      
    dlon *= 2.0/3.0
    ddep *= 2.0/3.0
    
end # end loop over iterations

# output final locations
clat1 = flatbest1
clon1 = flonbest1
cdep1 = fdepbest1
clat2 = flatbest2
clon2 = flonbest2
cdep2 = fdepbest2
resol = (dlat/(2.0/3.0))*degkm

# compute distance between cluster centroids
cdist = sqrt((cdep2-cdep1)^2 + (map_distance(clat1,clon1,clat2,clon2))^2)
                
# compute residual between observed and predicted travel time
sdist1 = map_distance(clat1.+qlat1,clon1.+qlon1,slat,slon)
sdist2 = map_distance(clat2.+qlat2,clon2.+qlon2,slat,slon)
@inbounds for ii=1:npick
    tt1 = ttTABs[iph[ii]](sdist1[ii],cdep1+qdep1[ii])
    tt2 = ttTABs[iph[ii]](sdist2[ii],cdep2+qdep2[ii])
    pdif = (tt2 + qorg2[ii]) - (tt1 + qorg1[ii])
    resid[ii] = tdif[ii] - pdif - torgdif
end

# rms residual and median absolute residual
rms = sqrt.(mean(resid.^2))
rmed = mad(resid,normalize=false) # default is to equate to sdev

# return results
return clat1, clon1, cdep1, clat2, clon2, cdep2, cdist, torgdif, resid, rms, rmed, resol

end

##### difclust2: L2 residual norm
function difclust2(qlat0::Float64,qlon0::Float64,qdep0::Float64,
    tdif::Vector{Float64},iph::Vector{Int8},
    slat::Vector{Float64},slon::Vector{Float64},
    qlat1::Vector{Float64},qlon1::Vector{Float64},
    qdep1::Vector{Float64},qorg1::Vector{Float64},
    qlat2::Vector{Float64},qlon2::Vector{Float64},
    qdep2::Vector{Float64},qorg2::Vector{Float64},
    ttTABs,boxwid::Float64,nit::Int64,degkm::Float64)

# initialize some variables
flatbest1, flonbest1, fdepbest1 = 0.0, 0.0, 0.0
flatbest2, flonbest2, fdepbest2 = 0.0, 0.0, 0.0
torgdif = 0.0

# extract npick
npick = length(tdif)
resid = zeros(npick)

# initialize box
dlat0 = 0.0
dlon0 = 0.0
ddep0 = 0.0
dlat = 0.5*boxwid/degkm
cosqlat = cosd(qlat0)
dlon = dlat/cosqlat

# careful near surface (added 07/2018)
zboxwid = boxwid
if (qdep0 < 0.0)
    zboxwid = min(1.0, boxwid)
elseif (qdep0 < boxwid/2.0)
    zboxwid = max(1.0,2.0*qdep0)
end
ddep = 0.5*zboxwid

# start iterations
for it=1:nit
    
    # grid search over box for best fit
    fitbest = 1.0e20
    tbest = 0.0
    
    # get trial locations: f1, f2 (3x3 grid for this box iteration)
    for iy = -1.0:1.0
        flat1 = qlat0 + dlat0 + dlat*iy
        flat2 = qlat0 - dlat0 - dlat*iy
        for ix = -1.0:1.0
            flon1 = qlon0 + dlon0 + dlon*ix
            flon2 = qlon0 - dlon0 - dlon*ix
            for iz = -1.0:1.0
                fdep1 = qdep0 + ddep0 + ddep*iz
                fdep2 = qdep0 - ddep0 - ddep*iz
                
                # compute predicted travel time and residuals w/observed
                sdist1 = map_distance(qlat1.+flat1,qlon1.+flon1,slat,slon)
                sdist2 = map_distance(qlat2.+flat2,qlon2.+flon2,slat,slon)                
                @inbounds for ii=1:npick
                    tt1 = ttTABs[iph[ii]](sdist1[ii],qdep1[ii]+fdep1)
                    tt2 = ttTABs[iph[ii]](sdist2[ii],qdep2[ii]+fdep2)
                    pdif = (tt2 + qorg2[ii]) - (tt1 + qorg1[ii])
                    resid[ii] = tdif[ii] - pdif # accounts for otime adjustment
                end
                
                # compute fit with L2 norm
                residval = mean(resid)
                fit = sum((resid.-residval).^2)
                
                # update best fit
                if fit < fitbest
                    fitbest = fit
                    flatbest1 = flat1
                    flonbest1 = flon1
                    fdepbest1 = fdep1
                    flatbest2 = flat2
                    flonbest2 = flon2
                    fdepbest2 = fdep2
                    tbest = residval
                end
                
            end
        end
    end # end of grid search
    
    # update best position
    dlat0 = flatbest1 - qlat0
    dlon0 = flonbest1 - qlon0
    ddep0 = fdepbest1 - qdep0
    torgdif = tbest # otime difference: residual between observed and predicted tdif

    # shrink box by 2/3 each iteration
    dlat *= 2.0/3.0      
    dlon *= 2.0/3.0
    ddep *= 2.0/3.0
    
end # end loop over iterations

# output final locations
clat1 = flatbest1
clon1 = flonbest1
cdep1 = fdepbest1
clat2 = flatbest2
clon2 = flonbest2
cdep2 = fdepbest2
resol = (dlat/(2.0/3.0))*degkm

# compute distance between cluster centroids
cdist = sqrt((cdep2-cdep1)^2 + (map_distance(clat1,clon1,clat2,clon2))^2)
                
# compute residual between observed and predicted travel time
sdist1 = map_distance(clat1.+qlat1,clon1.+qlon1,slat,slon)
sdist2 = map_distance(clat2.+qlat2,clon2.+qlon2,slat,slon)
@inbounds for ii=1:npick
    tt1 = ttTABs[iph[ii]](sdist1[ii],cdep1+qdep1[ii])
    tt2 = ttTABs[iph[ii]](sdist2[ii],cdep2+qdep2[ii])
    pdif = (tt2 + qorg2[ii]) - (tt1 + qorg1[ii])
    resid[ii] = tdif[ii] - pdif - torgdif
end

# rms residual and median absolute residual
rms = sqrt.(mean(resid.^2))
rmed = mad(resid,normalize=false) # default is to equate to sdev

# return results
return clat1, clon1, cdep1, clat2, clon2, cdep2, cdist, torgdif, resid, rms, rmed, resol

end

##### difclust3: L3 (robus mean) residual norm
function difclust3(qlat0::Float64,qlon0::Float64,qdep0::Float64,
    tdif::Vector{Float64},iph::Vector{Int8},
    slat::Vector{Float64},slon::Vector{Float64},
    qlat1::Vector{Float64},qlon1::Vector{Float64},
    qdep1::Vector{Float64},qorg1::Vector{Float64},
    qlat2::Vector{Float64},qlon2::Vector{Float64},
    qdep2::Vector{Float64},qorg2::Vector{Float64},
    ttTABs,boxwid::Float64,nit::Int64,inorm::Int64,degkm::Float64)

# initialize some variables
flatbest1, flonbest1, fdepbest1 = 0.0, 0.0, 0.0
flatbest2, flonbest2, fdepbest2 = 0.0, 0.0, 0.0
torgdif = 0.0

# extract npick
npick = length(tdif)
resid = zeros(npick)

# initialize box
dlat0 = 0.0
dlon0 = 0.0
ddep0 = 0.0
dlat = 0.5*boxwid/degkm
cosqlat = cosd(qlat0)
dlon = dlat/cosqlat

# careful near surface (added 07/2018)
zboxwid = boxwid
if (qdep0 < 0.0)
    zboxwid = min(1.0, boxwid)
elseif (qdep0 < boxwid/2.0)
    zboxwid = max(1.0,2.0*qdep0)
end
ddep = 0.5*zboxwid

# start iterations
for it=1:nit
    
    # grid search over box for best fit
    fitbest = 1.0e20
    tbest = 0.0
    
    # get trial locations: f1, f2 (3x3 grid for this box iteration)
    for iy = -1.0:1.0
        flat1 = qlat0 + dlat0 + dlat*iy
        flat2 = qlat0 - dlat0 - dlat*iy
        for ix = -1.0:1.0
            flon1 = qlon0 + dlon0 + dlon*ix
            flon2 = qlon0 - dlon0 - dlon*ix
            for iz = -1.0:1.0
                fdep1 = qdep0 + ddep0 + ddep*iz
                fdep2 = qdep0 - ddep0 - ddep*iz
                
                # compute predicted travel time and residuals w/observed
                sdist1 = map_distance(qlat1.+flat1,qlon1.+flon1,slat,slon)
                sdist2 = map_distance(qlat2.+flat2,qlon2.+flon2,slat,slon)                
                @inbounds for ii=1:npick
                    tt1 = ttTABs[iph[ii]](sdist1[ii],qdep1[ii]+fdep1)
                    tt2 = ttTABs[iph[ii]](sdist2[ii],qdep2[ii]+fdep2)
                    pdif = (tt2 + qorg2[ii]) - (tt1 + qorg1[ii])
                    resid[ii] = tdif[ii] - pdif # accounts for otime adjustment
                end
                
                # compute fit (depends on norm)
                residval, fit = robomean(resid,0.1,10)
                
                # update best fit
                if fit < fitbest
                    fitbest = fit
                    flatbest1 = flat1
                    flonbest1 = flon1
                    fdepbest1 = fdep1
                    flatbest2 = flat2
                    flonbest2 = flon2
                    fdepbest2 = fdep2
                    tbest = residval
                end
                
            end
        end
    end # end of grid search
    
    # update best position
    dlat0 = flatbest1 - qlat0
    dlon0 = flonbest1 - qlon0
    ddep0 = fdepbest1 - qdep0
    torgdif = tbest # otime difference: residual between observed and predicted tdif

    # shrink box by 2/3 each iteration
    dlat *= 2.0/3.0      
    dlon *= 2.0/3.0
    ddep *= 2.0/3.0
    
end # end loop over iterations

# output final locations
clat1 = flatbest1
clon1 = flonbest1
cdep1 = fdepbest1
clat2 = flatbest2
clon2 = flonbest2
cdep2 = fdepbest2
resol = (dlat/(2.0/3.0))*degkm

# compute distance between cluster centroids
cdist = sqrt((cdep2-cdep1)^2 + (map_distance(clat1,clon1,clat2,clon2))^2)
                
# compute residual between observed and predicted travel time
sdist1 = map_distance(clat1.+qlat1,clon1.+qlon1,slat,slon)
sdist2 = map_distance(clat2.+qlat2,clon2.+qlon2,slat,slon)
@inbounds for ii=1:npick
    tt1 = ttTABs[iph[ii]](sdist1[ii],cdep1+qdep1[ii])
    tt2 = ttTABs[iph[ii]](sdist2[ii],cdep2+qdep2[ii])
    pdif = (tt2 + qorg2[ii]) - (tt1 + qorg1[ii])
    resid[ii] = tdif[ii] - pdif - torgdif
end

# rms residual and median absolute residual
rms = sqrt.(mean(resid.^2))
rmed = mad(resid,normalize=false) # default is to equate to sdev

# return results
return clat1, clon1, cdep1, clat2, clon2, cdep2, cdist, torgdif, resid, rms, rmed, resol

end
