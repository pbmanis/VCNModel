TITLE Ohmic Leak Current
:
: Ohmic leak current that can be set on a per-section basis
:
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
		(S) = (siemens)
}

NEURON {
	SUFFIX leak
	NONSPECIFIC_CURRENT i
	RANGE i, erev, gbar
	}
PARAMETER {
	gbar = 2e-5	(S/cm2)  < 0, 1e9 >
	erev = -50	(mV)
}
ASSIGNED {
	i  (mA/cm2)
	v  (mV)
}
BREAKPOINT {
	i = gbar*(v - erev)
}
