# constants.py: A python file containing useful constants for Project Sparkles
# v.1: 19/11/2021 - Alejandro S. Borlaff - a.s.borlaff@nasa.gov

# Useful definitions 
# Photometry: 
# http://faraday.uwyo.edu/~admyers/ASTR5160/handouts/516016.pdf
# https://www.mso.anu.edu.au/~amedling/obstech/obstech_5_fluxes.pdf
# ESA Land Observations: https://earth.esa.int/landtraining09/D1Lb1_PotuckovaOpticalBasics.pdf
# https://arxiv.org/pdf/1507.03578.pdf

# Moon and Earth parameters:
# https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html

# Physical Constants 
R_earth = 6371E+3 # Radius of the Earth, m
R_moon = 1738.2E+3 # Radius of the Moon, m (Allen 1976)
d_moon_earth = 384400E+3 # Distance Moon - Earth, m 
h_sat = 600E+3 # Altitude of the spacecraft, m
cs = 299792458. # Speed of light, m/s
h = 6.62607004E-34 # Planck's constant m2 kg / s
albedo_earth = 0.434 # Average value for the albedo of the Earth 
albedo_moon = 0.12 # Average value for the albedo of the Moon
d_sun_earth = 149597870E+3 # m

# Statistical constants 
sigma1=0.682689492137086 
sigma2=0.954499736103642 
sigma3=0.997300203936740 
