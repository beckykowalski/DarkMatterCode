import numpy as np

##### define mathematical converstion factors

# 1 GeV/c^2 = 1.7827 * 10^-27 kg
GeVtokg = 1.78266192 * 10.**(-27.)

# 1 keV = 1.602 * 10^-16 J
keVtoJ = 1.602176 * 10.**(-16.)

# 1 GeV = 1.602 * 10^-10 J
GeVtoJ = 1.602176 * 10.**(-10.)

MeVtoJ = 1.602176 * 10.**(-13.)

# 1 km = 1000 m
kmtom = 1000.

# 1km = 100000 cm
kmtocm = 100000.

# 1 square meter = 10000 square cm
m2tocm2 = 10000.

# second to day = 60 sec in min * 60 min in hour * 24 hours
sectoday = 1./ (60.*60.*24.)
sectoyear = 1./ (60.*60.*24.*365)

# velocities 
vE = 230. # km / s, Earth's nominal velocity
vEsc = 544. # km / s, Galactic Escape veloicty
v0 = 238. # km / s, Maxwellian DM velocity distribution
vEavg = 29.8 # km / s , Earth's average galactocentric speed 

k0 = (np.pi * v0*v0)**(3./2.)
rho = 0.3 # GeV / (c^2 cm^3), DM local density

# Earth's orbit vectors, taken from https://arxiv.org/pdf/2105.00599.pdf
timeomega = 0.0172 # units: days^-1
earthT0 = 0.218 # unit years: should be March 22, 2018. Can make unit of days

EpsilonSpringEquinox = [-0.0504, 0.4946, -0.8677] # March 21
EpsilonSummerSolstice = [-0.9941, -0.1088, -0.0042] # June 21

# Solar peculiar velocities, taken from https://arxiv.org/pdf/2105.00599.pdf
## units in km/s
vSolar_R = 11.1
vSolar_Theta = 12.2
vSolar_Phi = 7.3

## numerical variables
pi = 3.1415926
twopi = 6.2831853
pi2 = 9.8696044
sqrtpi = 1.772454
sqrt_32 = 1.22474
sqrt_23 = 0.816496581
sqrt_2pi_over1 = 0.398942
OneThird = 0.333333
OneSixth = 0.166667
Five_18ths = 0.277778



