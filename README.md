# AGN-X-ray-Variability

xabsgrid_ufo: Absorption spectra generated by the xabs model (https://spex-xray.github.io/spex-help/models/xabs.html).

ioa_functions: Functions needed for the interpolation routine.

F_var_functions: functions needed to generate the excess variance spectra.

ioa_bil_int_all: Bilinear interpolation code. Input x,y (logxi,lognh) and recieve flux at this energy and ionisation/column density, using flux at nearest (logxi,lognh) points in the .fits file for the interpolation. Genearates new flux at (logxi,lognh) and energy, which is then iterated all energies to produce the RMS spectra.

F_var_ie_interval: For constant column density, takes 10 values of flux at values of ie energy in rangeie_l to ie_h. [we find these values of flux at energy ind_e (index of energy in .fits dataset) and (ie,col) using interpolation from ioa_functions]. Calc F_var for each energy bin from this.[note: energy range has been cropped between indexes 4422 and 6554 (0.5 to 10keV) as we are interested in the iron lines at ~7keV].

F_var_col_interval: Same routine as F_var_ie_interval but for constant ionisation energy and a range of columnn densities.

F_var_ie_maximise: Calculates RMS spectra for n ranges of ionisation energies between imax and imin, with the ranges having intervals of step. Ratio measurement calculates the relative strength of the Fe line and plots the figure with parameters which maximise its amplitude. Ratios and ionisation ranges are also returned. Column density is constant.

F_var_col_maximise: Same routine as F_var_ie_maximise but with ranges of column densities and constant ionisation. 


[Note that the maximise routines may plot graphs which have high absorption activity at lower energies instead of those which maximise the iron line, these outputs can be disregarded and the next highest ratio used instead].
