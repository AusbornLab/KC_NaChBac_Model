: NaChBac channel model based on Strege et al., 2023.
: Model developed by: Anthony Moreno-Sanchez 


: Set units
UNITS {
   (S)  = (siemens)
   (mV) = (millivolt)
   (mA) = (milliamp)
}

: Declare the channel
: the suffix declaration in the below section refers to the name that this mod file is called by in python/hoc code (see line 20 of the python script)
NEURON {
   SUFFIX na_bac_strege
   USEION na READ ena WRITE ina
   RANGE gmax, ina, m , h, ina2
}

: Set default parameters
: this is where you adjust gmax
PARAMETER {
   gmax = 0.12 (S/cm2) : the value to the left is the value in Neuron's original hh.mod file
}

ASSIGNED {
   ena (mV)
   ina (mA/cm2)
   v  (mV)
   ina2 ()
}

STATE { m h }

: Calculate activation and current
: this is where you adjust degrees for activation/inactivation
BREAKPOINT {
   SOLVE states METHOD cnexp
   ina = gmax * m^3 * h * (v - ena)
   ina2 =  m^3 * h 
}

: Initialize activation to steady state
: this is where steady state curves are defined
INITIAL { 

   : values extracted from Strege et al. are tuned by cmaes
   m = 1.0 / (1.0 + exp( (-47.20712608-v) / 8.171491024))
   h = 1.0 / (1.0 + exp( (-56.68326979-v) / -6.078198844 ))
}

: Calculate change in activation 
: equations should have same values as ones above, make sure that the derivative states you are using matches the steady-state equations you are using 
: this means that the values in your m' should match the values in your m, and same with your h' and h
DERIVATIVE states {
   UNITSOFF

   m' = (1.0 / (1.0 + exp( (-47.20712608-v) / 8.171491024 )) - m) / taum(v)
   h' = (1.0 / (1.0 + exp( (-56.68326979-v) / -6.078198844 )) - h) / tauh(v)

   UNITSON
}

: Function for activation time constant
: taum is defined here
FUNCTION taum(v (mV)) (/ms) {
   UNITSOFF

   taum = -125253.86 / (1 + exp(-0.06 * (v - -187.78))) + 125258.35


   UNITSON
}

: Function for inactivation time constant
: tauh is defined here
FUNCTION tauh(v (mV)) (/ms) {
   UNITSOFF
   
   tauh = 36.11944425 / ( 1 + exp((-24.58946004-v)/-5.15842944) ) + 165.144429082

   UNITSON
}
