# LapseRates

Makes plots of lapse rates, moist adiabats, and climate sensitivity. Change the parameters at the beginning of the file to set temperature and pressure ranges, decide which plots are shown and whether to save them, etc.

Moist adiabatic lapse rates are given in K/Km and calculated using the Clausius-Clapeyron relation and hydrostatic balance. Note: here I assume the total pressure is equal to the pressure of non-condensible gas - appropriate for the dilute limit where the saturation mixing ratio is much less than one.

Moist adiabatic temperature profiles are calculated by integrating the moist adiabatic lapse rate.

Climate sensitivity is calculated here as change in surface temperature for a change in radiative forcing, assuming outgoing radiation is emmitted from one level with a characteristic emission pressure and temperature. It does not include greenhouse forcing from water vapor / the condensible species in comparing the moist and dry cases.

This project was originally created for MIT 12.315 Atmospheric Radiation and Convection
