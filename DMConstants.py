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
vEsc = 544. # km / s, Galactic Escape veloicty
v0 = 238. # km / s, Maxwellian DM velocity distribution
vEavg = 29.8 # km / s , Earth's average galactocentric speed 

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

# Earth's orbit vectors, taken from https://arxiv.org/pdf/2105.00599.pdf

epsilon1_R = -0.0505
epsilon1_Theta = 0.4944
epsilon1_Phi = -0.8677
epsilon2_R = -0.9941
epsilon2_Theta = -0.1088
epsilon2_Phi = -0.0042
timeometa = 0.0172 # units: days^-1
deltaT = 0.218 # unit years: should be March 22, 2018. Can make unit of days

EpsilonSpringEquinox = [-0.0505, 0.4944, -0.8677] # March 21
EpsilonSummerSolstice = [0.9441, 0.1088, 0.0042] # June 21

# Solar peculiar velocities, taken from https://arxiv.org/pdf/2105.00599.pdf
## units in km/s

vSolar_R = 11.1
vSolar_Theta = 12.2
vSolar_Phi = 7.3




