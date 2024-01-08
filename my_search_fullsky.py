import pandas as pd
import numpy as np
from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table, vstack, hstack, join
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
from shared import save_spectrum

sdss_data_release = 16

def search_sdss_fullsky():
    query = f'''
    select
        p.objid, e.specObjID, p.type, e.ra, e.dec, e.z as redshift, p.u,p.g,p.r,p.i,p.z, 
        e.V_Ha_6562, e.Flux_Ha_6562, e.EW_Ha_6562, e.Amplitude_Ha_6562,
        e.V_Hb_4861, e.Flux_Hb_4861, e.EW_Hb_4861, e.Amplitude_Hb_4861,
        e.V_OIII_4958, e.Flux_OIII_4958, e.EW_OIII_4958, e.Amplitude_OIII_4958,
        e.V_OIII_5006, e.Flux_OIII_5006, e.EW_OIII_5006, e.Amplitude_OIII_5006,
        e.V_SII_6716, e.Flux_SII_6716, e.EW_SII_6716, e.Amplitude_SII_6716,
        e.V_SII_6730, e.Flux_SII_6730, e.EW_SII_6730, e.Amplitude_SII_6730
    FROM PhotoObj AS p
    JOIN emissionLinesPort as e on e.specObjID = p.specObjID
    LEFT JOIN ROSAT ON p.objid = ROSAT.OBJID
    WHERE 
        e.V_Ha_6562 > 500 and e.Flux_Ha_6562 > 1000 and
        e.V_Hb_4861 > 500 and e.Flux_Hb_4861 > 100 and
        e.V_OIII_4958 > 500 and e.Flux_OIII_4958 > 100 and
        e.V_OIII_5006 > 500 and e.Flux_OIII_5006 > 300 and
        e.V_SII_6716 > 500 and e.Flux_SII_6716 > 50 and
        e.V_SII_6730 > 500 and e.Flux_SII_6730 > 30 and
        ROSAT.OBJID IS NULL
    '''
    xid = SDSS.query_sql(query, data_release=sdss_data_release)
    return xid

result = search_sdss_fullsky()
rows = len(result)
print(f'SDSS search returns {rows} objects')

def compute_median_flux(fits_data, z):
    hdu = fits_data[0]
    # Access the data part of the FITS file
    data = hdu[1].data
    # The spectral flux is usually in a field named 'flux' and the wavelength in 'loglam' or 'wavelength'
    flux = data['flux']
    wavelength = 10**data['loglam']  # If loglam is given, convert to linear scale
    median_flux_regions = [4850 * (1+z), 5050 * (1+z)]
    region_mask = (wavelength >= median_flux_regions[0]) & (wavelength <= median_flux_regions[1])
    median_flux = np.median(flux[region_mask])
    return median_flux

# append a new column to indicate the signal
result['HII_Singal'] = 0

for i in tqdm(range(rows - 1, -1, -1)):
    row = result[i]
    specObjID = row['specObjID']
    z = row['redshift']
    query = f'SELECT TOP 1 specObjID, run2d, plate, fiberID, mjd FROM SpecObj WHERE specObjID={specObjID}'
    result2 = SDSS.query_sql(query, data_release=sdss_data_release)
    should_exclude = True
    if result is not None:
        try:
            fits_data = SDSS.get_spectra(matches=result2, data_release=sdss_data_release)
            if fits_data is not None:
                median_flux = compute_median_flux(fits_data, z)
                signal_flux = (row['Flux_Hb_4861'] + row['Flux_OIII_4958'] + row['Flux_OIII_5006']) / 3
                if signal_flux / median_flux > 15:
                    row['HII_Singal'] = signal_flux / median_flux
                    should_exclude = False
        except:
            pass
    if should_exclude:
        result.remove_row(i)

print(f'After removing the low signal rows, the table has {len(result)} objects')
result.write('data/my_search_result_fullsky.csv', format='csv', overwrite=True)
print('CSV saved')
