: Transient sodium channel

: Set units
UNITS {
   (S)  = (siemens)
   (mV) = (millivolt)
   (mA) = (milliamp)
}

: Declare transient sodium channel
NEURON {
   SUFFIX nat
   USEION na READ ena WRITE ina
   RANGE natgmax, ina, m , h, taum, tauh, ina2nat
}

: Set default parameters
PARAMETER {
   natgmax = 0.0012 (S/cm2)
}

ASSIGNED {
   ina (mA/cm2)
   ena (mV)
   v (mV)
   ina2nat ()
}

STATE { m h }

: Calculate current for transient sodium channel
BREAKPOINT {
   SOLVE states METHOD cnexp
   ina = natgmax * m^3 * h * (v - ena)
   ina2nat = m^3 * h
}

: Initialize activation and inactivation to steady state
INITIAL { 
m = 1.0 / (1.0 + exp( 0.1121*(-29.13- v)))
h = 1.0 / (1.0 + exp( -0.2*(-47- v)))
}

: Calculate change in activation / inactivation 
DERIVATIVE states {
   UNITSOFF
m' = ((1.0 / (1.0 + exp( 0.1121*(-29.13- v)))) - m) / taum(v)
h' = ((1.0 / (1.0 + exp( -0.2*(-47- v)))) - h) / tauh(v)

   UNITSON
}

: Function for activation time constant
FUNCTION taum(Vm (mV)) (/ms) {
   UNITSOFF
   :taum = -3.43 / (1 + exp(-0.17 * (v - (-46.03)))) + 3.56
   taum = 0.1270 + 3.434 / (1 + exp((v + 45.35) / 5.98))
   UNITSON
}

: Function for inactivation time constant
FUNCTION tauh(Vm (mV)) (/ms) {
   UNITSOFF
   tauh = 0.36 + exp( (v + 20.65) / (-10.47) )
   UNITSON
}