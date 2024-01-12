from astropy.io import fits

flare_fits = 'GOES_XRS_historical_finalversion.fits'

fitsfile = fits.open(flare_fits)
data = fitsfile[1].data
header = fitsfile[1].header

print(data.columns)
		