# toolbox.py: A python file containing useful functions and constants for Project Sparkles
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
s1_down_q = (1-sigma1)/2
s1_up_q = 1 - s1_down_q
s2_down_q = (1-sigma2)/2
s2_up_q = 1 - s2_down_q
s3_down_q = (1-sigma3)/2
s3_up_q = 1 - s3_down_q

############################################################ 
############### Photometry Functions ####################### 
############################################################ 

# Cheatsheet
# AB mag in W/Hz/m2: mAB = -2.5*log10(f) - 56.1  

def flux2mag(flux): 
    # Input flux in W/m2/Hz
    # Output magAB 
    # Ref: http://faraday.uwyo.edu/~admyers/ASTR5160/handouts/516016.pdf
    return(-2.5*np.log10(flux) - 56.1)

def mag2flux(mag): 
    # Output flux in W/m2/Hz 
    # Input magAB 
    # Ref: http://faraday.uwyo.edu/~admyers/ASTR5160/handouts/516016.pdf
    return(10**((-56.1-mag)/2.5))

def flux_freq2wave(flux, lambda_ref):
    # Fλ  = Fν c / λ2  
    # Input: flux in W/m2/Hz, lambda_ref in m
    # Output: flux in W/m2/m
    return(cs*flux/lambda_ref**2)

def flux_wave2freq(flux, lambda_ref):
    # Fν    =    Fλ    λ2/c   
    # Input: flux in W/m2/m, lambda_ref in m
    # Output: flux in W/m2/Hz
    return(flux*lambda_ref**2/cs)

def angular_radius_earth_from_orbit(h_sat):
    # Input - Orbit altitude
    # Output - Apparent angular size of the Earth from orbit in degrees
    return(np.degrees(np.arctan(R_earth/(R_earth+h_sat))))


def sb_sunlit_moon(mag_sun, d_moon_obs):
    # Input: 
    # mag_sun: Apparent magnitude of the sun  
    # d_moon_observer: Distance of the observer to the surface of the Earth
    
    # This is the flux that arrives to Earth, the Solar Constant 
    flux_sun_at_moon = mag2flux(mag_sun) # W / m2 / Hz - We can aprox that it is the same for Moon and Earth
    
    # The Earth cross section is R_earth**2 * pi
    # Then the light has to illuminate half of the sphere until R_moon + h, times the albedo
    sunflux_reflected_by_moon = albedo_moon*flux_sun_at_moon*R_moon**2/(2*(R_moon+d_moon_obs)**2)
    
    # The angular radius of the earth in degrees as defined by the function 
    angular_radius_moon_from_obs = np.degrees(np.arctan(R_moon/(R_moon+d_moon_obs)))
    
    # The area is just r**2 * pi, times 3600**2 to switch from deg**2 to arcsec**2
    area_moon_from_obs = np.pi*angular_radius_moon_from_obs**2*3600**2
    
    # The surface brightness intensity is:
    sb_intensity = sunflux_reflected_by_moon/area_moon_from_obs
    return(sb_intensity)


def sb_sunlit_earth(mag_sun, d_earth_obs):
    # Input: 
    # mag_sun: Apparent magnitude of the sun  
    # albedo_earth: Albedo of the Earth
    # d_earth_observer: Distance of the observer to the surface of the Earth
    
    # This is the flux that arrives to Earth, the Solar Constant 
    flux_sun_at_earth = mag2flux(mag_sun) # W / m2 / Hz
    
    # The Earth cross section is R_earth**2 * pi
    # Then the light has to illuminate half of the sphere until R_earth + h, times the albedo
    sunflux_reflected_by_earth = albedo_earth*flux_sun_at_earth*R_earth**2/(2*(R_earth+d_earth_obs)**2)
    
    # The angular radius of the earth in degrees as defined by the function 
    angular_radius_earth_from_obs = angular_radius_earth_from_orbit(d_earth_obs)
    
    # The area is just r**2 * pi, times 3600**2 to switch from deg**2 to arcsec**2
    area_earth_from_obs = np.pi*angular_radius_earth_from_obs**2*3600**2
    
    # The surface brightness intensity is:
    sb_intensity = sunflux_reflected_by_earth/area_earth_from_obs
    return(sb_intensity)


def sb_moonlit_earth(mag_sun, d_earth_obs):
    # Input: 
    # mag_sun: Apparent magnitude of the sun  
    # d_earth_observer: Distance of the observer to the surface of the Earth
    
    # This is the flux that arrives to Earth, the Solar Constant 
    flux_sun_at_moon = mag2flux(mag_sun) # W / m2 / Hz
    
    # The Earth cross section is R_earth**2 * pi
    # Then the light has to illuminate half of the sphere until R_earth + h, times the albedo
    sunflux_reflected_by_moon_to_earth = albedo_moon*flux_sun_at_moon*R_moon**2/(2*(d_moon_earth)**2)

    # The intensity received by Earth from the full moon is
    I_moon_earth = sunflux_reflected_by_moon_to_earth*np.pi*R_earth**2

    # The flux of the night side of the Earth illuminated by the Moon at an altitude h_sat is:
    f_earthshine_moonlit = albedo_earth*I_moon_earth/(2*np.pi*(R_earth + d_earth_obs)**2)    
    
    # The angular radius of the earth in degrees as defined by the function 
    angular_radius_earth_from_obs = angular_radius_earth_from_orbit(d_earth_obs)
    
    # The area is just r**2 * pi, times 3600**2 to switch from deg**2 to arcsec**2
    area_earth_from_obs = np.pi*angular_radius_earth_from_obs**2*3600**2
    
    # The surface brightness intensity is:
    sb_intensity = f_earthshine_moonlit/area_earth_from_obs
    return(sb_intensity)
    
    
    
############################################################ 
############### Generic    Functions ####################### 
############################################################ 

def save_fits(array, name, header=None):
    hdu = fits.PrimaryHDU(header=header, data=array)
    hdul = fits.HDUList([hdu])
    os.system("rm " + name)
    hdul.writeto(name)
