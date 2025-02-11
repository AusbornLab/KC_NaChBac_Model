: Persistent sodium channel

: Set units
UNITS {
   (S)  = (siemens)
   (mV) = (millivolt)
   (mA) = (milliamp)
}

: Declare Persistent sodium channel
NEURON {
   SUFFIX nap
   USEION na READ ena WRITE ina
   RANGE napgmax, ina, m , h, taum, tauh, ina2nap
}

: Set default parameters
PARAMETER {
   napgmax = 0.0012 (S/cm2)
}

ASSIGNED {
   ina (mA/cm2)
   ena (mV)
   v (mV)
   ina2nap ()
}

STATE { m h}

: Calculate current for Persistent sodium channel
BREAKPOINT {
   SOLVE states METHOD cnexp
   ina = napgmax * m * (v - ena)
   ina2nap = m
}

: Initialize activation and inactivation to steady state
INITIAL { 
m = 1.0 / (1.0 + exp(0.2717 * (-48.77-v)))
}

: Calculate change in activation / inactivation 
DERIVATIVE states {
   UNITSOFF
   m' = (1.0 / (1.0 + exp(0.2717 * (-48.77-v))) - m) / taum(v)
   UNITSON
}

FUNCTION taum(Vm (mV)) (/ms) {
   taum = 1
}

