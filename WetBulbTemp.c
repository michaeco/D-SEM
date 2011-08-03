/*
 * SUMMARY:      WetBulbTemp.c - Calculate the wet bulb temperature
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Matthew Wiley
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:              mwwiley@u.washington.edu
 * ORIG-DATE:     11/03/2004
 * DESCRIPTION:  Calculates the wet bulb temperature given the relative
 *               humidity, air temperature, and pressure.  The wet bulb
 *               temperature is the temperature an air parcel would have 
 *               if cooled adiabatically to saturation at constant pressure
 *               by evaporation of water into it, all latent heat being 
 *               supplied by the parcel (1).  The temperature of an 
 *               evaporating cloud droplet, or a raindrop is equal to the 
 *               wet bulb temperature (2).  The raindrop temperature is used 
 *               in the stream temperature / heat balance calculations.  It
 *               may also potentially be used for a improved rain versus 
 *               snow algorithm.
 *
 * (1) Huschke, R. E., Ed., 1959: Glossary of Meteorology. Amer. Meteor. Soc., 638 pp.
 * (2) Handbook of Hydrology, David Maidment editor, p3.4
 *
 *
 * FUNCTIONS:    CalcWetBulbTemp()
 *               CalcSVP()
 * COMMENTS:
 *     
 */

#include <stdlib.h>
#include <math.h>
#include "settings.h"
/*********************************************************************************/
/* Calculates the Saturated Vapor Pressure , returns value in millibars
 *       References: Shuttleworth, W.J., Evaporation,  In: Maidment, D. R. (ed.),
 *                   Handbook of hydrology,  1993, McGraw-Hill, New York, etc..
 *                   Bras, R. A., Hydrology, an introduction to hydrologic
 *                   science, Addisson Wesley, Inc., Reading, etc., 1990.
 */

float CalcSVP(float T)
{
  float Pressure;
  Pressure = 610.78 * exp((double) ((17.269 * T) / (237.3 + T)))/100;
  return Pressure;
}


/**********************************************************************************/
/* Calculate Dew Point from supplied Relative Humidity, and Vapor Pressure */

float CalcDewPoint(float RH, float es)
{
  float Dp,e;
  e = es * RH/100;
  Dp = (243.5 * log(e/6.112))/(17.67 - log(e/6.112));
  return Dp;
}

/***********************************************************************************/
 /* Iribarne, J. V., and W. L. Godson, 1981: Atmospheric Thermodynamics. 3d ed. D. Reidel, 259 pp. */

float CalcWetBulbTemp(float RH, float Tair, float Press)
{

	float Td, Tw;  /* DewPoint and Wetbulb Temperatures */
	float deriv;
	float WetBulb = Tair;
	float e,  diff,   es , ew, ed, s;//peq,tcur,vpcur,
	int k = 0;
	float C15 = 26.66082;
	float C1 = 0.0091379024;
	float C2 = 6106.396;
	float F = 0.0006355; /* Cp/L*epsilon) (1/K) */

	Press /= 100;  /* convert from PA to millibar */
  	e = CalcSVP(Tair);
	Td = CalcDewPoint(RH,e);
	Tair += 273.15;
	Td += 273.15;

    es = exp(C15 - C1 * (Tair) - C2 / (Tair));
    ed = exp(C15 - (C1 * Td) - (C2 / Td));
    s = (es - ed) / (Tair - Td);
	Tw = ( Tair * F * Press + Td * s ) / (F * Press + s);
	ew = exp(C15 - (C1 * Tw) - (C2 / Tw));
	diff  = F * Press * (Tair - Tw ) - (ew - ed);
	while (diff > 0.00001 && k < 1000 ) {
                k++;
		ew = exp(C15 - (C1 * Tw) - (C2 / Tw));
		diff  = F * Press * (Tair - Tw ) - (ew - ed);
		deriv = ew * (C1 - (C2/(Tw*Tw))) - F * Press;
		WetBulb = Tw - 273.15;
		Tw = Tw - diff/deriv;
	}
    
	return WetBulb;
}

/**  For testing algorithm, not part of DHSVM   **/

/* 
int main() {
	float k,RH, Tair, Press, WetBulb;
	printf("Enter RH Tair Press: ");
        k = scanf("%f %f %f",&RH, &Tair, &Press);
	WetBulb = CalcWetBulbTemp(RH,Tair,Press);
        printf("WetBulb = %f\n", WetBulb);
	return 1;
}
*/

