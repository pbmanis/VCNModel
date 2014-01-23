TITLE CaPCalyx.mod  The presynaptic calcium current at the MNTB calyx of Held

COMMENT

NEURON implementation of Calcium current from Borst and Sakmann, 1998.
Equations are basic HH; parameters are taken from Figure 8 legend in that paper.

Original implementation by Paul B. Manis, October 2007

Contact: pmanis@med.unc.edu

ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
	SUFFIX CaPCalyx
	USEION ca READ eca WRITE ica
	RANGE gcapbar, gcap, ica
	GLOBAL minf, taum, alpha, beta
}


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 22 (degC)  : model is defined on measurements made at 23-24 deg in Germany
        dt (ms)
        eca = 43.9 (mV)
        gcapbar = 0.01 (mho/cm2) <0,1e9> : target is 48.9 nS total
		alpha (1/ms)
		beta (1/ms)
}

STATE {
        m
}

ASSIGNED {
 	ica (mA/cm2)
    gcap (mho/cm2)
    minf (1)
	taum (ms) 
	}

LOCAL mexp

BREAKPOINT {
	SOLVE states
    
	gcap = gcapbar*m*m
    ica = gcap*(v - eca)

}

UNITSOFF

INITIAL {
    trates(v)
    m = minf
}

PROCEDURE states() {  :Computes state variables m
	trates(v)      :             at the current v and dt.
	m = m + mexp*(minf-m)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

	q10 = 3^((celsius - 24)/10) : if you don't like room temp, it can be changed!

    minf = (1 / (1 + exp(-(v + 23.2) / 9.1)))

    alpha = 1.78*exp(v/23.3)
    beta = 0.140*exp(-v/15.0)
    taum = 1/(alpha + beta)
}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL minc
	TABLE minf, mexp
	DEPEND dt, celsius FROM -150 TO 150 WITH 300

    rates(v)    : not consistently executed from here if usetable_hh == 1
        : so don't expect the tau values to be tracking along with
        : the inf values in hoc

	minc = -dt * q10
	mexp = 1 - exp(minc/taum)
	}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON
