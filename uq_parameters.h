#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED

// Reversal potentials (mV)
// #define vL -62.5
// #define vNa 45.0
// #define vK -105.0
// #define vH -35.0
// #define vCa 120.0

// Reversal potentials (mV)
#define vL -62.5  // 17.5 -9.5 -35.5 -56.5 -68.5
#define vNa 45.0
#define vK -105.0
#define vH -35.0
#define vCa 120.0

// Conductances (mS/cm^2)
// #define gL 2.5
// #define gNa 29.17
// #define gK 12.96
// #define gH 20.0
// #define gCaT 11.0213 or 15.0213
// #define gNaP 8.3244
// #define gHVA 2.0
// #define gBK 5.0
// #define gMystery 10.0

#define gL 2.5
#define gNa 29.17
#define gK 12.96
// #define gH 20.0
// #define gCaT 15.0213
// #define gNaP 8.3244
#define gHVA 2.0
// #define gBK 5.0
#define gMystery 10.0

// Na
#define theta_mNa -25.0
#define sigma_mNa -6.5

// K
#define theta_nK -26.0
#define sigma_nK -9.0
#define tau_nK 10.0

// CaT
#define theta_mCaT -37.1
#define sigma_mCaT -4.8916
#define tau_mCaT 40.0
#define theta_hCaT -59.2
#define sigma_hCaT 11.2326
#define tau_hCaT 350.0

// NaP
#define theta_mNaP -40.0
#define sigma_mNaP -4.0
#define theta_hNaP -54.0
#define sigma_hNaP 5.0
#define tau_hNaP 500.0

// H
#define theta_mH -61.32
#define sigma_mH 5.855
#define tau_mH_T 100.0
#define delta_mH_T 0.205
#define theta_mH_T -65.95
#define sigma_mH_T 4.44

// BK
#define wBK_base 170 //base time constant (ms)

// Ca
#define Ca0 0.00002 // concentration of calcium (mM)
#define tau_Ca 8.0 //calcium concentration time constant (ms)
#define Ca_buffer 0.5 //accounts for quick calcium buffering
#define Ca_z 2.0 //unitless number
#define d 1.0 //depth where calcuim concentration is relevent (microns)

// other parameters
#define C 21.0 // capacitance (microF/cm^2)
#define F 96485.0 //Faraday's constant (C/mol)
// #define Iex -15 //-15 14 67 137 205

#endif
