"""
set up initial conditions for a MYStIX cluster, from real stellar positions and masses
guesses for velocities & velocity dispersions using GAIA data

"""
from __future__ import print_function

import os

import time

import numpy as np
import numpy.random as rnd

from numpy import pi, sqrt

from amuse import datamodel

from amuse.datamodel import AbstractParticleSet

from amuse.units import *

from amuse.io import write_set_to_file
from amuse.io import read_set_from_file

from amuse.ic.brokenimf import new_broken_power_law_mass_distribution
from amuse.ic.salpeter import new_salpeter_mass_distribution

from amuse.support.console import set_printing_strategy 

from amuse.ext.spherical_model import new_uniform_spherical_particle_distribution 

try:
    from usagi.ic.gasplummer import new_plummer_gas_model
    from usagi.ic.plummer import new_plummer_model
except:
    from gasplummer import new_plummer_gas_model
    from plummer import new_plummer_model

try:
    from usagi.plotting.starcluster import gas_stars_plot
except:
    from plotting import gas_stars_plot

from parameters import *
from argumentparser import new_IC_argument_parser


def new_xy_for_velocity(p):
    number_of_selected_items    = 0
    selected_values_for_x       = np.zeros(0)
    selected_values_for_y       = np.zeros(0)
    while ( number_of_selected_items < p.n_for_velocities ):
        x       = rnd.uniform(
                0.,
                1.0, 
                ( p.n_for_velocities - number_of_selected_items ),
                )
        y       = rnd.uniform(
                0.,
                0.1, 
                ( p.n_for_velocities - number_of_selected_items ),
                )
        g       = ( x**2 ) * np.power(
                1.0 - x**2,
                3.5,
                )
        compare = y <= g

        selected_values_for_x       = np.concatenate(
                (
                    selected_values_for_x, 
                    x.compress(compare),
                    ) 
                )
        selected_values_for_y       = np.concatenate(
                (
                    selected_values_for_x, 
                    y.compress(compare),
                    )
                )
        number_of_selected_items    = len( selected_values_for_x )
    return selected_values_for_x, selected_values_for_y

def new_velocities_spherical_coordinates(
        p, 
        radius,
        ):
    pi2         = pi * 2
    x, y        = new_xy_for_velocity(p)
    velocity    = x * sqrt(2.0) * np.power(
            1.0 + radius * radius,
            -0.25,)
    theta       = np.arccos(
            rnd.uniform(
                -1.0,
                1.0, 
                p.n_for_velocities,
                )
            )
    phi         = rnd.uniform(
            0.0,
            pi2, 
            p.n_for_velocities,
            )
    return(velocity, theta, phi)
 
def coordinates_from_spherical(
        p, 
        radius, 
        theta, 
        phi,
        ):
    x = radius * np.sin( theta ) * np.cos( phi )
    y = radius * np.sin( theta ) * np.sin( phi )
    z = radius * np.cos( theta )
    return (x,y,z)

def velocity_dispersion(particles):
    N                   = len(particles)
    dv                  = (
            particles.velocity - 
            particles.center_of_mass_velocity()
            )
    squarevelocities    = dv * dv
    sigma               = (
            (
                squarevelocities.sum() /
                ( N - 1 )
                ).sqrt()
            )
    return sigma

#################################################################################

class RealClusterIC (object):

    def __init__(
            self,
            p,
            ):
        if not p.seed is None:
            rnd.seed(p.seed)
 
        self.virial_ratio   = p.stars_virial_ratio

        self.initialize_stars()
        self.initialize_gas()

    # FIXME I am not happy about this thing being here...
    AbstractParticleSet.add_global_function_attribute(
            "velocity_dispersion",
            velocity_dispersion,
            )
#####################################################################      
    def initialize_stars(self):
        (ra,dec,vx,vy) = np.loadtxt(
                p.cluster_stars_file,usecols=(0,1,2,3),
                skiprows    = 1,
                unpack      = True,
                )
	
 
        nobsstars   = len(ra)
 
        rastars     = ra * pi / 180. #convert to radians
        decstars    = dec * pi / 180.
 
        (
                racentre, 
                deccentre, 
                rcmajor, 
                rcminor,
		paellipse,
                distance, 
                nstarstot,
                )   = np.genfromtxt(
                        p.cluster_file,
                        skip_header = p.cluster_number,
                        skip_footer = (
                            p.cluster_ntotal - 
                            p.cluster_number
                            ),
                        unpack      = True,
                        )
	
	#rcmajor and rcminor are in arcminutes, and are the axes of an ellipse 4 times the size of the fitted core. Want these to actually be the axes of the core, in parsec. Divide by 4, then convert to parsec using the cluster distance. 1 arcmin = 60 arcsec = 60/206265 radians

        rcmajor *= 0.25*60/206265*distance
	rcminor *= 0.25*60/206265*distance
	p.stars_n = nstarstot
	self.dist = distance

        # core radius is the harmonic mean of rcmajor and rc minor;
        # rscale = virial radius = 16 sqrt(2)/3 pi *rcore for a
        # Plummer sphere.

        p.rscale         = (
                2.0 / 
                (
                    1.0 / rcmajor + 
                    1.0 / rcminor
                    ) * 
                16. * sqrt(2) / 3. / pi 
                ) | units.parsec 
        p.variable_refresh()

        # Make converter, based on 1 MSun stars. It's ok.
        self.converter = nbody_system.nbody_to_si(
                p.stars_n | units.MSun, 
                p.rscale,
                )

        racentre        = racentre * pi / 180.  #convert to radians
        deccentre       = deccentre * pi / 180. 
        self.paellipse  = paellipse * pi / 180. 
 
        self.ellipticity    = ( rcmajor - rcminor ) / rcmajor
        self.yellip         = p.yellip*self.ellipticity
        self.xellip         = p.xellip*self.ellipticity
        self.zellip         = p.zellip*self.ellipticity
	 
        dcluster            = distance | units.parsec
 
        self.obsstars       = datamodel.Particles(nobsstars)
 
        imf_masses          = new_salpeter_mass_distribution(
                nobsstars,
                mass_max        = p.cluster_observed_mass_max,
                mass_min         = p.cluster_observed_mass_min,
                alpha          = -2.3
                )
        self.stellar_mass   = imf_masses.sum()
        self.obsstars.mass  = imf_masses
 
        dra     = rastars - racentre
        cosdec  = np.cos(deccentre)
        sindec  = np.sin(deccentre)
 
        x       = (
                np.cos(decstars) * np.sin(dra) / 
                (
                    cosdec * np.cos(decstars) * np.cos(dra) + 
                    sindec * np.sin(decstars)
                    )
                )
        y       = (
                (
                    cosdec * np.sin(decstars) - 
                    sindec * np.cos(decstars) * np.cos(dra)
                    ) /
                (
                    sindec * np.sin(decstars) + 
                    cosdec * np.cos(decstars) * np.cos(dra)
                    )
                )
      
        # Spherical trig -- everything is an angle, and everything is
        # in radians
        self.obsstars.x         = x * dcluster
        self.obsstars.y         = y * dcluster
        self.obsstars.z         = rnd.uniform(-p.obsstars_zfactor * rcminor,p.obsstars_zfactor  * rcminor,nobsstars,) | units.parsec

        self.obsstars.vx        = vx| units.kms 
        self.obsstars.vy        = vy| units.kms
        self.obsstars.vz        = rnd.uniform(-1.0,1.0,nobsstars,) | units.kms
        self.obsstars.radius    = p.stars_interaction_radius

#####################################################################
	
        # now do the background of lower-mass stars     
        nmorestars      = int(p.stars_n - nobsstars)
        imf_masses2     = new_broken_power_law_mass_distribution(
                nmorestars, 
                mass_boundaries = p.stars_mass_boundaries, 
                mass_max        = p.stars_mass_max, 
                alphas          = p.stars_mass_alphas, 
                random          = True,
                )
        othermass       = imf_masses2.sum()
        self.otherstars = new_plummer_model(
                nmorestars,
                xellip=self.xellip,
                yellip=self.yellip,
                zellip=self.zellip,
                convert_nbody   = self.converter,
                )
        self.otherstars.mass    = imf_masses2
        self.otherstars.move_to_center()
 
        # rotate the ellipse so that it lines up with the stars.
        # paellipse is measured from north (+y) axis towards the east
        # (+x axis)

        xnew    = (
                self.otherstars.x * np.cos(self.paellipse) + 
                self.otherstars.y * np.sin(self.paellipse)
                )
        ynew    = (
                self.otherstars.y * np.cos(self.paellipse) - 
                self.otherstars.x * np.sin(self.paellipse)
                )
 
	self.otherstars.x   = xnew
        self.otherstars.y   = ynew

	vxx=rnd.normal(-0.5488128,13.341549,nmorestars)
	vyy=rnd.normal(-2.285615,10.244671,nmorestars)
	
	vx1=vxx*4.74/distance
	vy1=vyy*4.74/distance

	self.otherstars.vx  = vx1 |units.kms
	self.otherstars.vy  = vy1 |units.kms
	self.otherstars.vz  = rnd.uniform(-1.0,1.0,nmorestars,) |units.kms
	
	self.otherstars.radius = p.stars_interaction_radius

########################################################################## 
        # put it all together
        self.stellar_mass  += othermass
 
        self.stars              = datamodel.Particles()
	self.stars.add_particles(self.obsstars)
        self.stars.add_particles(self.otherstars)
        self.stars.move_to_center()
        self.stars.scale_to_standard(
                convert_nbody   = self.converter,
               )

#########################################################################

    def initialize_gas(self):
        gas_mass    = self.stellar_mass * p.gas_fraction

        p.n_for_velocities     = int(gas_mass / p.gas_particle_mass)
 
        self.gas    = new_plummer_gas_model(
                p.n_for_velocities,
                xellip=self.xellip,
                yellip=self.yellip,
                zellip=self.zellip,
                convert_nbody   = self.converter,
                )
        self.gas.h_smooth   = self.converter.to_si(
                p.gas_smoothing_fraction
                )
        self.gas.u          = p.gas_u
        self.gas.mass       = p.gas_particle_mass 
        self.gas.move_to_center()
 
        # rotate the ellipse so that it lines up with the stars.
        # paellipse is measured from north (+y) axis towards the east
        # (+x axis)
        xnew    = (
                self.gas.x * np.cos(self.paellipse) + 
                self.gas.y * np.sin(self.paellipse)
                )
        ynew    = (
                self.gas.y * np.cos(self.paellipse) - 
                self.gas.x * np.sin(self.paellipse)
                )
 
        self.gas.x  = xnew
        self.gas.y  = ynew
 	
	vxx=rnd.normal(-0.5488128,13.341549,p.n_for_velocities)
	vyy=rnd.normal(-2.285615,10.244671,p.n_for_velocities)
	
	vx1=vxx*4.74/self.dist
	vy1=vyy*4.74/self.dist
 
 	self.gas.vx = vx1 |units.kms
	self.gas.vy = vy1 |units.kms
	self.gas.vz = rnd.uniform(-1.0,1.0,p.n_for_velocities)|units.kms

######################################################################################
 
    def write_data(self):
        write_set_to_file(
                self.stars,
                p.dir_initialconditions + p.stars_initial_file,
                'amuse',
                )
        write_set_to_file(
                self.gas,
                p.dir_initialconditions + p.gas_initial_file,
                'amuse',
                )
 
########################################################################################TOTAL CLUSTER

class AllClusterIC (object):

    def __init__(
            self,
            p,
            stars_sets,
            gas_sets,):
        self.stars_sets = stars_sets
        self.gas_sets   = gas_sets

        (
                racentre,
                deccentre,
                rcmajor,
                rcminor,
		paellipse,
                distance,
                nstarstot,
                )   = np.loadtxt(
                        p.cluster_file,
                        skiprows    = 1,
                        unpack      = True,
                        )
        rcmajor		*= 0.25*60/206265*distance
	rcminor		*= 0.25*60/206265*distance
	self.racentre   = racentre * pi / 180.
        self.deccentre  = deccentre * pi / 180.
       	self.dcluster   = distance[0] | units.parsec

########################################################

    def read_combine_clusters(
            self,
            ):
        self.allstars   = datamodel.Particles()
        self.allgas     = datamodel.Particles()
	
        cosdec  = np.cos(p.cluster_dec_avg)
        sindec  = np.sin(p.cluster_dec_avg)
        dra     = self.racentre - p.cluster_ra_avg
        dx      = (
                np.cos(self.deccentre) * np.sin(dra) /
                ( 
                    cosdec * np.cos(self.deccentre) * np.cos(dra) + 
                    sindec * np.sin(self.deccentre)
                    )
                )
        dy      = (
                (
                    cosdec * np.sin(self.deccentre) - 
                    np.cos(self.deccentre) * sindec * np.cos(dra)
                    ) /
                (
                    sindec * np.sin(self.deccentre) + 
                    cosdec * np.cos(self.deccentre) * np.cos(dra)
                    )
                )
        dx     *= self.dcluster
        dy     *= self.dcluster

# include a random variation in z as well. p.cluster_zmax gives the extent of the z variation, such that the centre of each cluster will be chosen to be between -1.0*p.cluster_zmax and p.cluster_zmax (in parsec)

        dz                   = rnd.uniform(
            -1.0*p.cluster_zmax.value_in(units.parsec),
            p.cluster_zmax.value_in(units.parsec),
            len(self.stars_sets),
            ) | units.parsec
        

        self.gas_mass   = 0.0 | units.MSun

        for i in range(0,len(self.stars_sets)):
           self.stars       = self.stars_sets[i]
           self.stars.x    += dx[i]
           self.stars.y    += dy[i]
           self.stars.z    += dz[i]
           self.allstars.add_particles(self.stars)

           self.gas         = self.gas_sets[i]
           self.gas.x      += dx[i]
           self.gas.y      += dy[i]
           self.gas.z      += dz[i]
           if i == 0:
              self.gasu = self.gas.u[0]
              self.gash = self.gas.h_smooth[0]
              print(self.gasu, self.gash)

           self.gas_mass   += self.gas.mass.sum()
           print(i, self.gas_mass)
           self.allgas.add_particles(self.gas)
                                                                       
        self.nclustered = len(self.allstars.x)
        self.allstars.move_to_center()
        self.allgas.move_to_center()

######################################################################  UNCLUSTERED

    def add_other_stars(self):
        
	(ra,dec,vx,vy) = np.loadtxt(
                p.cluster_stars_u,usecols=(0,1,2,3),
                skiprows    = 1,
                unpack      = True,
                )
        rastars             = ra * pi / 180.
        decstars            = dec * pi / 180.
        nstars              = len(ra)
	nunclustered        = int(nstars)
        unclustered_stars   = datamodel.Particles(nunclustered)
        imf_masses          = new_salpeter_mass_distribution(
                nunclustered, 
                mass_max        = p.cluster_observed_mass_max,
                mass_min         = p.cluster_observed_mass_min,
                alpha          = -2.3
                )
        unclustered_stars.mass		= imf_masses
	unclustered_stellar_mass	= imf_masses.sum()

        dra     = rastars - p.cluster_ra_avg
        cosdec  = np.cos(p.cluster_dec_avg)
        sindec  = np.sin(p.cluster_dec_avg)

	unclustered_converter = nbody_system.nbody_to_si(
			unclustered_stellar_mass,
			p.rscale,
			)

        x = (
                np.cos(decstars) * np.sin(dra) / 
                ( 
                    cosdec * np.cos(decstars) * np.cos(dra) + 
                    sindec * np.sin(decstars)
                    )
                )
        y = (
                (
                    cosdec * np.sin(decstars) - 
                    sindec * np.cos(decstars) * np.cos(dra)
                    ) / 
                (
                    sindec * np.sin(decstars) + 
                    cosdec * np.cos(decstars) * np.cos(dra)
                    )
                )
	
        unclustered_stars.x = x * self.dcluster 
	unclustered_stars.y = y * self.dcluster 
	unclustered_stars.z = rnd.uniform(
                -1.0*p.unclustered_zdistance.value_in(units.parsec),
                p.unclustered_zdistance.value_in(units.parsec),
                nstars,
                ) | units.parsec
	
	unclustered_stars.vx = vx|units.kms
	unclustered_stars.vy = vy|units.kms
	unclustered_stars.vz = rnd.uniform(-1.0,1.0,nstars,) | units.kms
        unclustered_stars.radius    = p.stars_interaction_radius

        self.allstars.add_particles(unclustered_stars)

############################################################----UNKNOWN STARS---

        (ra,dec,vx,vy)= np.loadtxt(
                p.cluster_stars_x,usecols=(0,1,2,3),
                skiprows    = 1,
                unpack      = True,
                )
        rastars             = ra * pi / 180.
        decstars            = dec * pi / 180.
        nstars              = len(ra)
        nunknown	    = int(nstars)
        unknown_stars       = datamodel.Particles(nunknown)
        imf_masses          = new_salpeter_mass_distribution(
                nunknown, 
                mass_max        = p.cluster_observed_mass_max,
                mass_min         = p.cluster_observed_mass_min,
                alpha          = -2.3
                )
        unknown_stars.mass  = imf_masses
	unknown_stellar_mass = imf_masses.sum()

        dra     = rastars - p.cluster_ra_avg
        cosdec  = np.cos(p.cluster_dec_avg)
        sindec  = np.sin(p.cluster_dec_avg)

        unknown_converter = nbody_system.nbody_to_si(
			unknown_stellar_mass,
			p.rscale,
			)

	x   = (
                np.cos(decstars) * np.sin(dra) /
                (
                    cosdec * np.cos(decstars) * np.cos(dra) +
                    sindec * np.sin(decstars)
                    )
                )
        y   = (
                (
                    cosdec * np.sin(decstars) - 
                    sindec * np.cos(decstars) * np.cos(dra)
                    ) / 
                (
                    sindec * np.sin(decstars) + 
                    cosdec * np.cos(decstars) * np.cos(dra)
                    )
                )
	
        unknown_stars.x = x * self.dcluster
	unknown_stars.y = y * self.dcluster
	unknown_stars.z = rnd.uniform(
                -1.0*p.unclustered_zdistance.value_in(units.parsec),
                p.unclustered_zdistance.value_in(units.parsec), 
		nunknown,
		) |units.parsec

	unknown_stars.vx = vx |units.kms
	unknown_stars.vy = vy |units.kms
	unknown_stars.vz = rnd.uniform(-1.0,1.0,nstars,) | units.kms

        unknown_stars.radius    = p.stars_interaction_radius

        self.allstars.add_particles(unknown_stars)

#################################################################LOW MASS STARS

# Number of low-mass unclustered/unknown stars calculated from extending the IMF (assuming broken power law with alpha = 1.3 between 0.2 and 0.5 Msun; and alpha = 2.3 between 0.5 and 100 Msun)
# notes from 13 June 2016
	
	nlow = int(nunclustered * (5.838 - p.cluster_observed_mass_min.value_in(units.MSun)**(-1.3))/ (p.cluster_observed_mass_min.value_in(units.MSun)**(-1.3) - p.cluster_observed_mass_max.value_in(units.MSun)**(-1.3)))
	lowmass_stars   = datamodel.Particles(nlow)
        imf_masses      = new_broken_power_law_mass_distribution(
                nlow, 
                mass_boundaries = p.stars_mass_boundaries,
                mass_max        = p.stars_mass_max, 
                alphas          = p.stars_mass_alphas, 
                random          = True,
                )
        othermass         = imf_masses.sum()
	
	lowmass_converter = nbody_system.nbody_to_si(
			othermass,
			p.rscale,
			)

	lowmass_stars   = new_uniform_spherical_particle_distribution(
                nlow, 
                p.lowmass_radius, 
                othermass, 
                type    = "random",
                )

        lowmass_stars.mass      = imf_masses

	vxx=rnd.normal(-0.5488128,13.341549,nlow,)
	vyy=rnd.normal(-2.285615,10.244671,nlow,)
	
	d=self.dcluster.value_in(units.parsec)
	vx1=vxx*4.74/d
	vy1=vyy*4.74/d

	lowmass_stars.vx  = vx1 |units.kms
	lowmass_stars.vy  = vy1 |units.kms
	lowmass_stars.vz  = rnd.uniform(-1.0,1.0,nlow,) |units.kms
	 

        lowmass_stars.radius    = p.stars_interaction_radius
        self.allstars.add_particles(lowmass_stars)

############################################################FILAMENT

    def add_filament(self):

	p.gas_filament_n    = int(p.gas_filament_mass / p.gas_particle_mass)

        gas_converter   = nbody_system.nbody_to_si(
                self.gas_mass,
                p.filament_rscale,
                )

        filament    = new_plummer_gas_model(
                p.gas_filament_n,
                xellip=p.filament_xellip, 
		yellip=p.filament_yellip,
                zellip=p.filament_zellip, 
		convert_nbody   = gas_converter
                )
        filament.h_smooth   = self.gash
        filament.mass       = p.gas_particle_mass 
	filament.u          = self.gasu
	
        paellipse   = p.filament_paellipse * pi / 180.
        
	xnew        = (
                filament.x * np.cos(paellipse) + 
                filament.y * np.sin(paellipse)
                )
        ynew        = (
                filament.y * np.cos(paellipse) - 
                filament.x * np.sin(paellipse)
                )
        filament.x  = xnew
        filament.y  = ynew


	vxx=rnd.normal(-0.5488128,13.341549,p.gas_filament_n)
	vyy=rnd.normal(-2.285615,10.244671,p.gas_filament_n)
	
	d=self.dcluster.value_in(units.parsec)
	vx1=vxx*4.74/d
	vy1=vyy*4.74/d

	filament.vx  = vx1 |units.kms
	filament.vy  = vy1 |units.kms
	filament.vz  = rnd.uniform(-1.0,1.0,p.gas_filament_n,) |units.kms
	 
	 
	
        self.allgas.add_particles(filament)

##################################################################

    def write_all(self):
        write_set_to_file(
                self.allstars,
                p.dir_initialconditions + p.stars_initial_file,
                'amuse',
                )
        write_set_to_file(
                self.allgas,
                p.dir_initialconditions + p.gas_initial_file,
                'amuse',
                )

#########################################################################------MAIN PROGRAM--------

if __name__ in("__main__"):
    p       = Parameters()
    p       = new_IC_argument_parser(p)
    
    set_printing_strategy(
            "custom",
            preferred_units = p.log_units_preferred,
            precision       = p.log_units_precision,
            )

    timestamp   = time.strftime("%Y%m%d-%H%M%S")

    stars_sets   = []
    gas_sets    = []

    for n in range(1, p.cluster_ntotal+1):
        p.cluster_number        = n
        print("Creating initial conditions for %s cluster %s...  "%(
                p.cluster_datarelease,
                p.cluster_label[p.cluster_number-1],
                ), end=' ')
        p.cluster_stars_file    = p.dir_input + "/%s%sstars.txt"%(
                p.cluster_datarelease,
                p.cluster_label[p.cluster_number-1],
                )

        p.dir_simulation    = "./Results-%s/Run-%s/"%(
                timestamp,
                p.cluster_label[p.cluster_number-1],
                )

        p.variable_refresh()

        if not os.path.exists(p.dir_simulation):
            os.makedirs(p.dir_simulation)
            os.makedirs(p.dir_initialconditions)
            os.makedirs(p.dir_plots)

        system = RealClusterIC(
                p,
                )

        system.write_data()
        p.plot_axes_x = "x"
        p.plot_axes_y = "y"
        i = gas_stars_plot(
            0,
            0.0 | units.Myr,
            system.gas,
            system.stars,
            p
            )
        p.plot_axes_x = "x"
        p.plot_axes_y = "z"
        j = gas_stars_plot(
            0,
            0.0 | units.Myr,
            system.gas,
            system.stars,
            p
            )
        p.plot_axes_x = "y"
        p.plot_axes_y = "z"
        k = gas_stars_plot(
            0,
            0.0 | units.Myr,
            system.gas,
            system.stars,
            p
            )
        print("Done!")

        stars_sets.append(system.stars)
        gas_sets.append(system.gas)


    #### From here, we build the full set ###

    p.dir_simulation    = "./Results-%s/All/"%(
            timestamp,
            )

    p.variable_refresh()

    if not os.path.exists(p.dir_simulation):
        os.makedirs(p.dir_simulation)
        os.makedirs(p.dir_initialconditions)
        os.makedirs(p.dir_plots)

#######################################################################

    full_system = AllClusterIC(
            p,
            stars_sets,
            gas_sets,
    )

    full_system.read_combine_clusters()
    full_system.add_other_stars()
    full_system.add_filament()
    full_system.write_all()

    p.plot_minx = -4.0 | units.parsec
    p.plot_maxx =  4.0 | units.parsec
    p.plot_miny = -4.0 | units.parsec
    p.plot_maxy =  4.0 | units.parsec

    p.plot_axes_x = "x"
    p.plot_axes_y = "y"
    i = gas_stars_plot(
        0,
        0.0 | units.Myr,
        full_system.allgas,
        full_system.allstars,
        p
        )

    p.plot_axes_x = "x"
    p.plot_axes_y = "z"
    j = gas_stars_plot(
        0,
        0.0 | units.Myr,
        full_system.allgas,
        full_system.allstars,
        p
        )

    p.plot_axes_x = "y"
    p.plot_axes_y = "z"
    k = gas_stars_plot(
        0,
        0.0 | units.Myr,
        full_system.allgas,
        full_system.allstars,
        p
        )
