import numpy as np
import astropy
import healpy as hp

from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy import constants as cte
from astropy import units as u
from matplotlib.colors import LogNorm
from astropy import wcs
from astropy.time import Time



def index_coords(data, origin=None):

    """
    A function that creates a map of the numbered pixels from an origin selected.
    
    Parameters
    -------
    data: 2-d numpy array.
        The map that you want to number.
    origin: optional, [x0, y0].
        Origin selected.
    """
    
    ny, nx = data.shape[:2]
    if origin is None:
        origin_x, origin_y = nx // 2, ny // 2
    else:
        origin_y, origin_x = origin
        if origin_y < 0:
            origin_y += ny
        if origin_x < 0:
            origin_x += nx

    x, y = np.meshgrid(np.arange(float(nx)) - origin_x,
                       origin_y - np.arange(float(ny)))
    return x, y
    
    
def map_of_coordinates(coords, n_pix, plate_scale):

    """
    A function that calculates the coordinates of every pixel of a detector. It is returned a map of coordinates (ra, dec) in degrees (SkyCoord (ICRS)).
    
    Parameters
    -------
    coordinates: SkyCoord (ICRS): (ra, dec) in deg.
            Coordinates of the pointing of the telescope.
    plate_scale: float (in arcsec/pixel).
            Plate scale of the telescope.
    n_pix: float.
        Number of pixels in x-axis of the detector
    """
    
    
    # Create a new WCS object.  The number of axes must be set
    # from the start
    data = np.zeros([n_pix,n_pix])

    xcenter = int(n_pix/2)
    ycenter = int(n_pix/2)

    x,y = index_coords(data, origin=(xcenter, ycenter))

    cdelt=np.array([-1.,1.])/3600 * plate_scale
    crpix = np.array([0., 0.])

    w = wcs.WCS(naxis=2)
   
    pa = 0
    crpix = [0., 0.]
    w.wcs.crpix = crpix
    w.wcs.crval = [coords.ra.deg, coords.dec.deg]
    w.wcs.cdelt = cdelt
    w.wcs.crota = [0, -pa]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    rac_det, dec_det = w.all_pix2world(x, y, 1)

    list_detector = SkyCoord(ra=rac_det, dec=dec_det, unit='deg', frame='icrs')
    
    return list_detector




def earth_ndi_detector(ndi_function, sat_distance, plate_scale, n_pix):

    """
    A function that calculates the straylight from the earth in a LEO telescope.
    
    Parameters
    -------
    ndi_function: function.
            A function dependent on the angular separetion of the telescope.
    sat_distance: float (in km).
            Altitude of the satellite (orbit).
    plate_scale: float (in arcsec/pixel).
            Plate scale of the telescope.
    n_pix: float.
        Number of pixels in x-axis of the detector
            
    Note:
        The function returns a map in W/m^2/Hz units.
    """
    
    

    def indextoradec(index):
        theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
        return np.degrees(np.pi*2.-phi), -np.degrees(theta-np.pi/2.)

    def radectoindex(ra,dec):
        return hp.pixelfunc.ang2pix(NSIDE,np.radians(-dec+90.),np.radians(360.-ra))


    #################
    ### CONSTANTS
    #################

    R_earth = 6378.0 # km
    # Earth surface brightness
    Rt = R_earth *1e3
    h = sat_distance*1e+3 #  m
    m_sun = -26.74 # solar magnitude
    m_moon = -12.7 # moon magnitude
    albedo_earth = 0.3
    filter_width = 3000 # A

    earth_theta = np.arctan(Rt/(Rt+h))*2
    print(np.degrees(earth_theta))

    print("Area earth")
    area_earth = (np.degrees(earth_theta)/2)**2*np.pi # deg^2
    print(area_earth)


    sun_int = 10**(-0.4*(56.1+m_sun))   # W / m2 / Hz
    moon_int = 10**(-0.4*(56.1+m_moon)) # W / m2 / Hz

    print("Intensity received at the surface of Earth from Moon: W / m2 / Hz")
    print(moon_int)

    night_earth_int = moon_int * albedo_earth

    I_moon_dunes = night_earth_int*Rt**2/(4*(h + Rt)**2) # W / m2 / Hz
    
    print("Intensity received at the satellite from Earth illuminated by Moon: W / m2 / Hz")
    print(I_moon_dunes)

    mag_earth_illuminated_by_moon = -2.5*np.log10(I_moon_dunes) - 56.1
    mu_earth_illuminated_by_moon = -2.5*np.log10(I_moon_dunes/(area_earth*3600*3600)) - 56.1
    
    print("Surface brightness of the Earth illuminated by Moon: mags/arcsec2")
    print(mu_earth_illuminated_by_moon)

    ##################

    day_earth_int = (moon_int+sun_int) * albedo_earth

    I_day_dunes = day_earth_int*Rt**2/(4*(h + Rt)**2) # W / m2 / Hz
    
    print("Intensity received at the satellite from Earth illuminated by Sun: W / m2 / Hz")
    print(I_day_dunes)

    mag_earth_illuminated_by_daylight = -2.5*np.log10(I_day_dunes) - 56.1
    mu_earth_illuminated_by_daylight = -2.5*np.log10(I_day_dunes/(area_earth*3600*3600)) - 56.1

    print("Surface brightness of the Earth illuminated by Sun: mags/arcsec2")
    print(mu_earth_illuminated_by_daylight)
    
    
    # Coordinates of the center of the Earth in the POV of the satellite
    lon = 0
    lat = np.pi/2
    vec = hp.ang2vec(lat, lon)
    
    # Create a healpix sphere:
    NSIDE = 32
    NPIX = hp.nside2npix(NSIDE)
    
    # Create a circle from the size of the earth watched at the satellite:
    ipix_disc = hp.query_disc(NSIDE, vec=vec, radius=np.radians(np.arcsin(R_earth/(R_earth+sat_distance))*180/np.pi))

    # Divide the circle in day and night
    lon_list_day, lat_list_day = hp.pix2ang(nside=2**5, ipix=ipix_disc[int(len(ipix_disc)/2):,])
    lon_list_night, lat_list_night = hp.pix2ang(nside=2**5, ipix=ipix_disc[:int(len(ipix_disc)/2)])

    # RA and DEC of the pixels
    ra_night, dec_night = indextoradec(ipix_disc[:int(len(ipix_disc)/2)])
    ra_day, dec_day = indextoradec(ipix_disc[int(len(ipix_disc)/2):])
    coords_list_night = SkyCoord(ra_night, dec_night, frame='icrs', unit = 'deg')
    coords_list_day = SkyCoord(ra_day, dec_day, frame='icrs', unit = 'deg')
    
    # Coords of the pointing of the telescope in the healpix sphere:
    satellite_coords = SkyCoord(np.degrees(np.pi*2.-0), -np.degrees(0-np.pi/2.), frame='icrs', unit = 'deg')

    list_detector = map_of_coordinates(satellite_coords, n_pix, plate_scale):

    separation_night =  np.zeros([len(coords_list_night),n_pix,n_pix])
    separation_day =  np.zeros([len(coords_list_day),n_pix,n_pix])
    
    for i in range(len(coords_list_night)):
        separation_day[i]=coords_list_day[i].separation(list_detector).deg
        separation_night[i]=coords_list_night[i].separation(list_detector).deg


    surface_pixels_arcsec_2 = hp.nside2pixarea(NSIDE, degrees=True)*3600**2 # arcsecs^2 of each pixel
    
    mag_night_pixel = mu_earth_illuminated_by_moon - 2.5*np.log10(surface_pixels_arcsec_2)
    mag_day_pixel = mu_earth_illuminated_by_daylight - 2.5*np.log10(surface_pixels_arcsec_2)

    flux_night = 10**((mag_night_pixel + 56.1)/-2.5) # W/m2/Hz
    flux_day = 10**((mag_day_pixel + 56.1)/-2.5) # W/m2/Hz
    
    flux_detector_day = ndi_function(separation_day)*flux_day
    flux_detector_night = ndi_function(separation_night)*flux_night
    
    total_night = np.sum(flux_detector_night, axis=0)
    total_day = np.sum(flux_detector_day, axis=0)
    
    total = total_day + total_night
    
    return total

