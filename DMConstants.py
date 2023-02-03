import numpy as np

##### define mathematical converstion factors

# 1 GeV/c^2 = 1.7827 * 10^-27 kg
GeVtokg = 1.7827 * 10.**(-27.)

# 1 keV = 1.602 * 10^-16 J
keVtoJ = 1.602 * 10.**(-16.)

# 1 GeV = 1.602 * 10^-10 J
GeVtoJ = 1.602 * 10.**(-10.)

# 1 km = 1000 m
kmtom = 1000.

# 1km = 100000 cm
kmtocm = 100000.

# 1 square meter = 10000 square cm
m2tocm2 = 10000.

# second to day = 60 sec in min * 60 min in hour * 24 hours
sectoday = 1./ (60.*60.*24.)
sectoyear = 1./ (60.*60.*24.*365)

##### define constants

vE = 230. # km / s, Earth's nominal velocity
vEsc = 600. # km / s, Galactic Escape veloicty
v0 = 238. # km / s, Maxwellian DM velocity distribution
##### define mathematical converstion factors

# 1 GeV/c^2 = 1.7827 * 10^-27 kg
GeVtokg = 1.7827 * 10.**(-27.)

# 1 keV = 1.602 * 10^-16 J
keVtoJ = 1.602 * 10.**(-16.)

# 1 GeV = 1.602 * 10^-10 J
GeVtoJ = 1.602 * 10.**(-10.)

# 1 km = 1000 m
kmtom = 1000.

# 1km = 100000 cm
kmtocm = 100000.

# 1 square meter = 10000 square cm
m2tocm2 = 10000.

# second to day = 60 sec in min * 60 min in hour * 24 hours
sectoday = 1./ (60.*60.*24.)
sectoyear = 1./ (60.*60.*24.*365)

##### define constants

vE = 230. # km / s, Earth's nominal velocity
vEsc = 600. # km / s, Galactic Escape veloicty
#v0 = 220. # km / s, Maxwellian DM velocity distribution
v0 = 238. # km / s, Maxwellian DM velocity distribution
#### define mathematical converstion factors

# 1 GeV/c^2 = 1.7827 * 10^-27 kg
GeVtokg = 1.7827 * 10.**(-27.)

# 1 keV = 1.602 * 10^-16 J
keVtoJ = 1.602 * 10.**(-16.)

# 1 GeV = 1.602 * 10^-10 J
GeVtoJ = 1.602 * 10.**(-10.)

# 1 km = 1000 m
kmtom = 1000.

# 1km = 100000 cm
kmtocm = 100000.

# 1 square meter = 10000 square cm
m2tocm2 = 10000.

# second to day = 60 sec in min * 60 min in hour * 24 hours
sectoday = 1./ (60.*60.*24.)
sectoyear = 1./ (60.*60.*24.*365)

##### define constants

# for Maxwellian Velocity Distribution 
vE = 230. # km / s, Earth's nominal velocity
vEsc = 600. # km / s, Galactic Escape veloicty
v0 = 238. # km / s, Maxwellian DM velocity distribution
k0 = (np.pi * v0*v0)**(3./2.)


rho = 0.3 # GeV / (c^2 cm^3), DM local density
