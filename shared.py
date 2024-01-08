import matplotlib.pyplot as plt

def save_spectrum(sp, z, title, filename):
    # Assuming 'sp' is the FITS object returned by SDSS.get_spectra()
    # Open the FITS file (the first spectrum in case there are multiple)
    hdu = sp[0]

    # Access the data part of the FITS file
    data = hdu[1].data

    # The spectral flux is usually in a field named 'flux' and the wavelength in 'loglam' or 'wavelength'
    flux = data['flux']
    wavelength = 10**data['loglam']  # If loglam is given, convert to linear scale

    # H_alpha
    Ha_rest = 6562.8
    Ha_observed = Ha_rest * (1 + z)

    # H_beta
    Hb_rest = 4861
    Hb_observed = Hb_rest * (1 + z)

    # OIII
    O3_rest = 5006
    O3_observed = O3_rest * (1 + z)

    # SII
    S2_rest = 6730
    S2_observed = S2_rest * (1 + z)

    # Plot the spectrum
    plt.figure(figsize=(10, 6))
    plt.plot(wavelength, flux, drawstyle='steps-mid')
    plt.xlabel(f'Wavelength (Å), z = {z}')
    plt.ylabel('Flux (10^-17 erg/s/cm^2/Å)')
    plt.title(title)

    # Add a vertical line for the Hα wavelength
    plt.axvline(x=Ha_observed, color='r', linestyle='dotted', linewidth=1.0, label=f'Ha line (6563))')
    plt.axvline(x=Hb_observed, color='cyan', linestyle='dotted', linewidth=1.0, label=f'Hb line (4861)')
    plt.axvline(x=O3_observed, color='b', linestyle='dotted', linewidth=1.0, label=f'O III line (5006)')
    plt.axvline(x=S2_observed, color='orange', linestyle='dotted', linewidth=1.0, label=f'S II line (6730)')

    # Add a legend
    plt.legend()

    plt.grid(True)
    plt.savefig(filename)
    plt.close()