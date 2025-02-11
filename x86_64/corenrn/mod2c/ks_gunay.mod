: Slow potassium channel

: Set units
UNITS {
   (S)  = (siemens)
   (mV) = (millivolt)
   (mA) = (milliamp)
}

: Declare transient sodium channel
NEURON {
   SUFFIX ks_gunay
   USEION k READ ek WRITE ik
   RANGE ksgmax, ik, m , h, taum, tauh, ik2
}

: Set default parameters
PARAMETER {
   ksgmax = 0.0012 (S/cm2)
}

ASSIGNED {
   ik (mA/cm2)
   ek (mV)
   v (mV)
   ik2 ()
}

STATE { m h }

: Calculate current for transient sodium channel
BREAKPOINT {
   SOLVE states METHOD cnexp
   ik = ksgmax * m^4 * (v - ek)
   ik2 = m^4 
}

: Initialize activation and inactivation to steady state
INITIAL { 
m = 1.0 / (1.0 + exp(0.0502 * (-12.85-v)))
}

: Calculate change in activation / inactivation 
DERIVATIVE states {
   m' = ((1.0 / (1.0 + exp(0.0502 * (-12.85-v)))) - m) / taum(v)
}

: Function for activation time constant
FUNCTION taum(Vm (mV)) (/ms) {
   taum = 2.03 + 1.96 / (1 + exp( (Vm - 30.83) / 3.12))
}

