from parcels import FieldSet, ParticleSet, Variable, JITParticle, ErrorCode, AdvectionRK4, convert_IndexedOutputToArray
from parcels import random
from datetime import timedelta as delta
import datetime
import numpy as np
from glob import glob
import math


def set_hycom_grid(t=None):
    t0 = datetime.date(2016, 7, 1)
    ddir = '/Volumes/data02/HYCOMdata/GLBa0.08_expt90_surf/hycom_GLBu0.08_912_'
    if t is None:
        files = [ddir + (t0 - delta(days=i)).strftime('%Y%m%d') + "00_t000.nc" for i in range(2, -1, -1)]
    else:
        files = [ddir + (t0 - delta(days=t)).strftime('%Y%m%d') + "00_t000.nc"]
    print t, files
    filenames = {'U': files, 'V': files}
    variables = {'U': 'u', 'V': 'v'}
    dimensions = {'lat': 'Latitude', 'lon': 'Longitude', 'time': 'MT', 'depth': 'Depth'}
    indices = {'lat': range(1850, 2001), 'lon': range(3700, 4100)}
    fset = FieldSet.from_netcdf(filenames, variables, dimensions, indices=indices)
    return fset


def set_startlocs():
    locs = [[35.29311, 32.93944],
            [35.32631, 32.93527],
            [35.36705, 32.92333],
            [35.33255, 33.48277],
            [35.33463, 33.49305],
            [35.35416, 33.59750],
            [35.41191, 33.83416],
            [35.41592, 33.86361],
            [35.54833, 34.17166],
            [35.60072, 34.33388],
            [35.62633, 34.36972],
            [35.66666, 34.57222],
            [35.64116, 34.54694],
            [35.52297, 34.33972],
            [35.365, 34.14],
            [35.279, 33.98],
            [35.168,  33.97]]
    locs = map(list, zip(*locs))
    return locs[0], locs[1]


def OutOfBounds(particle, fieldset, time, dt):
    particle.delete()


def AddAge(particle, fieldset, time, dt):
    particle.age = particle.age + math.fabs(dt)
    if particle.age > fieldset.maxage:
        particle.delete()


def run_globcurrent_particles():
    fset = set_hycom_grid()
    fset.maxage = 365*86400
    i0 = len(fset.U.time)

    lats, lons = set_startlocs()
    lons = [lon + 360 for lon in lons]

    class PlasticParticle(JITParticle):
        age = Variable('age', dtype=np.int32, initial=0)

    pset = {}
    for i in range(len(lons)):
        pset[i] = ParticleSet.from_list(fieldset=fset, pclass=PlasticParticle, lon=[lons[i]],
                                        lat=[lats[i]], time=fset.U.time[-1])
        pset[i].pfile = pset[i].ParticleFile(name='plastic_cyprusbeaches_s'+"%0.3d"%i, type="indexed")
        pset[i].pfile.write(pset[i], pset[i][0].time)

    kernels = pset[0].Kernel(AdvectionRK4) + AddAge
    for t in range(i0, 730, 1):
        for pkey, p in pset.items():
            p.execute(kernels, starttime=p[0].time, runtime=delta(days=1), dt=-delta(minutes=5),
                      recovery={ErrorCode.ErrorOutOfBounds: OutOfBounds})
            if t<365:
                p.add(PlasticParticle(lon=lons[pkey], lat=lats[pkey], fieldset=fset))
            p.pfile.write(p, p[0].time)

        fset.advancetime(set_hycom_grid(t))

    for p in pset.values():
        convert_IndexedOutputToArray(file_in=p.pfile.name + ".nc", file_out=p.pfile.name + "_array.nc")


run_globcurrent_particles()
