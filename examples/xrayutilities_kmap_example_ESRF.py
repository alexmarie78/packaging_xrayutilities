# This file is part of xrayutilities.
#
# xrayutilities is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) 2014 Raphael Grifone <raphael.grifone@esrf.fr>
# Copyright (C) 2014 Dominik Kriegner <dominik.kriegner@gmail.com>

import os

import xrayutilities as xu
from matplotlib.pylab import *

import xrayutilities_id01_functions as id01

# define some convenience variables
en = 9000.0  # x-ray energy in eV
home = "\\homefolder"
workdir = os.path.join(home, 'work')
specdir = home  # location of spec file
sample = "sample"  # sample name -> used as spec file name
ccdfiletmp = "ccdfilename"
# region of interest on the detector; useful to reduce the amount of data
roi = [0, 516, 50, 300]

# define experimental geometry and detector parameters
# 2S+2D goniometer (simplified ID01 goniometer, sample mu,phi detector nu,del
qconv = xu.experiment.QConversion(['z+', 'y-'], ['z+', 'y-'], [1, 0, 0])
# convention for coordinate system: x downstream; z upwards; y to the
# "outside" (righthanded)
hxrd = xu.HXRD([1, 1, 0], [0, 0, 1], en=en, qconv=qconv)
hxrd.Ang2Q.init_area('z-', 'y+', cch1=200.07, cch2=297.75, Nch1=516, Nch2=516,
                     pwidth1=9.4489e-05, pwidth2=9.4452e-05, distance=1.,
                     detrot=-0.801, tiltazimuth=30.3, tilt=1.611, roi=roi)


def hotpixelkill(ccd):
    """
    function to remove hot pixels from CCD frames or apply any other filter.
    """
    # insert your own code here
    return ccd

U = numpy.identity(3)  # orientation matrix of the sample

scannr = numpy.arange(0, 100)

nx, ny = 100, 50

# specfile, scannumbers, nx,ny, motornames, optional column names (ID01
# values are default)
fss = xu.io.FastScanSeries(sample, scannr, nx, ny, 'Mu', 'Eta', 'Nu', 'Delta',
                           xmotor='adcY', ymotor='adcX', ccdnr='imgnr',
                           path=specdir)
# retrace clean all scans
fss.retrace_clean()
# align different scans (misalignment needs to be determined manually or if
# your application allows also by a 2d correllation code
# see e.g. scipy.signal.correlate2d
# fss.align(deltax,deltay)

# real space grid
g2d = xu.Gridder2D(nx, ny)
# read all motor positions from the data files
fss.read_motors()

# plot counter intensity on a grid
for idx, fs in enumerate(fss.fastscans):
    print(idx)
    g2d.Clear()
    g2d(fs.data['adcX'], fs.data['adcY'], fs.data['mpx4int'])
    plt.figure()
    plt.contourf(g2d.xaxis, g2d.yaxis, g2d.data.T, 50)
    plt.xlabel("X (um)")
    plt.ylabel("Y (um)")

posx, posy = 50, 30
# reduce data: number of pixels to average in each detector direction
nav = [2, 2]

qnx, qny, qnz = (80, 100, 101)
g = fss.gridRSM(posx, posy, qnx, qny, qnz, qconv, ccdfiletmp, path='',
                roi=roi, nav=nav, typ='index', filterfunc=hotpixelkill,
                UB=U)
# with typ='real' the position should be real space coordinates. with
# typ='index' the posx,y should specify indices within the range(nx)
# range(ny)

# g now contains a Gridder3D object which can be used for visualization
# see xrayutilities_example_plot_3D_ESRF_ID01 for example code

# Note: if you instead want to grid in two dimensions you can decrease one
# of qnx,qny or qnz to 1 and interprate the data in 2D
# Copyright (C) 2014,2018 Dominik Kriegner <dominik.kriegner@gmail.com>

# 3S+2D goniometer (ID01 goniometer, sample mu, eta, phi detector
# nu, del, mpxy, mpxz
# convention for coordinate system: x downstream; z upwards; y to the
# "outside" (righthanded)
# QConversion will set up the goniometer geometry. So the first argument
# describes the sample rotations, the second the detector rotations and the
# third the primary beam direction.
# For this consider the following right handed coordinate system (feel free to
# use your conventions):
# x: downstream (direction of primary beam)
# y: out of the ring
# z: upwards
# The outer most sample rotation (so the one mounted on the floor) is one
# which turns left-handed (-) around the z-direction -> z- (mu)
# The second sample rotation ('eta') is lefthanded (-) around y -> y-
qconv = xu.experiment.QConversion(['z-', 'y-', 'z-'],
                                  ['z-', 'y-', 'ty', 'tz'],
                                  [1, 0, 0])
hxrd = xu.HXRD([1, 1, 0], [0, 0, 1], qconv=qconv, sampleor='z+')
hxrd._A2QConversion.init_area('z-', 'y+', cch1=333.94, cch2=235.62, Nch1=516,
                              Nch2=516, pwidth1=5.5000e-02, pwidth2=5.5000e-02,
                              distance=0.53588*1000, detrot=-1.495,
                              tiltazimuth=155.0, tilt=0.745, Nav=(2, 2))

#############################################################
# load spec-file data
fnames = ('KMAP_2017_fast_00001.spec',)

scannrs = []
for fn in fnames:
    sf = xu.io.SPECFile(fn, path=datadir)
    scannrs.append(list(range(len(sf))))

nx, ny = 150, 150

# specfile, scannumbers, nx,ny, motornames, optional column names (ID01
# values are default)
fss = xu.io.FastScanSeries(fnames, scannrs, nx, ny, 'mu', 'eta', 'phi', 'nu',
                           'del', 'mpxy', 'mpxz', xmotor='adcY', ymotor='adcX',
                           ccdnr='imgnr', path=datadir)

#############################################################
# 3D RSM from summing over X,Y
# now all EDF files are parsed, this will take a while and some memory
qconv.energy = id01.getmono_energy(sf[0])
xu.config.VERBOSITY = 2
g3d = fss.get_average_RSM(81, 82, 83, qconv, datadir=id01.datadir,
                          replacedir=id01.repl_n, roi=None, nav=(4, 4),
                          filterfunc=deadpixelkill)
xu.config.VERBOSITY = 1
numpy.savez_compressed('RSM3D.npz', qx=g3d.xaxis, qy=g3d.yaxis, qz=g3d.zaxis,
                       data=g3d.data)


figure()
subplot(221)
pcolormesh(qx, qy, data.sum(axis=2).T, norm=mpl.colors.LogNorm())
xlabel(r'Qx ($\AA^{-1}$)')
ylabel(r'Qy ($\AA^{-1}$)')

subplot(222)
pcolormesh(qx, qz, data.sum(axis=1).T, norm=mpl.colors.LogNorm())
xlabel(r'Qx ($\AA^{-1}$)')
ylabel(r'Qz ($\AA^{-1}$)')

subplot(223)
pcolormesh(qy, qz, data.sum(axis=0).T, norm=mpl.colors.LogNorm())
xlabel(r'Qy ($\AA^{-1}$)')
ylabel(r'Qz ($\AA^{-1}$)')
tight_layout()

#############################################################
# 2D real space maps for selected region in Q-space
qr = [0.57, 0.62, -0.20, -0.16, 3.47, 3.50]  # [xmin, xmax, ymin, ..., zmax]
x, y, data = fss.get_sxrd_for_qrange(qr, qconv, datadir=id01.datadir,
                                     replacedir=id01.repl_n)
numpy.savez_compressed('output_sxrd_map.npz', x=x, y=y, data=data)


figure()
lev_exp = np.linspace(np.log10(data.min()),
                      np.log10(data.max()), 100)
levs = np.power(10, lev_exp)
tricontourf(y, x, data, levs, norm=mpl.colors.LogNorm())
axis('scaled')
ylabel('piy (um)')
xlabel('pix (um)')
tight_layout()
