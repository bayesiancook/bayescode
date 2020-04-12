#pragma once

// amino acids in one-letter code alphabetical order
// A C D E F G H I K L M N P Q R S T V W Y

// composition
const double grantham_com[] = {0,2.75,1.38,0.92,0,0.74,0.58,0,0.33,0,0,1.33,0.39,0.89,0.65,1.42,0.71,0,0.13,0.2};

// polarity
const double grantham_pol[] = {8.1,5.5,13.0,12.3,5.2,9.0,10.4,5.2,11.3,4.9,5.7,11.6,8.0,10.5,10.5,9.2,8.6,5.9,5.4,6.2};

// volume
const double grantham_vol[] = {31,55,54,83,132,3,96,111,119,111,105,56,32.5,85,124,32,61,84,170,136};

// associated weights
const double grantham_wcom = 1.833;
// 0.459 in Beaulieu et al
const double grantham_wcom_selac = 0.459;
const double grantham_wpol = 0.1018;
const double grantham_wvol = 0.000399;

// total distance between amino acids a and b: weighted sum of squared differences
/*
double tcom = grantham_com[b] - grantham_com[a];
double tpol = grantham_pol[b] - grantham_pol[a];
double tvol = grantham_vol[b] - grantham_vol[a];
double d = sqrt(grantham_wcom*tcom*tcom + grantham_wpol*tpol*tpol + grantham_wvol*tvol*tvol);
*/

