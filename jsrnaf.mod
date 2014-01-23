TITLE jsrnaf.mod VCN Na conductance, fast model

COMMENT
gnaf is the modified form used in his
1993 M.S. thesis (as in Rothman et al., J. Neurophysiol. 70:2562, 1993),
with rapid recovery from inactivation for potentials below rest.

Implementation by Paul B. Manis, April and Sept, 1999.
ENDCOMMENT

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(nA) = (nanoamp)
}

? interface
NEURON {
	SUFFIX na
	USEION na READ ena WRITE ina
	RANGE gnabar, gna, vsna 
	RANGE minf, hinf, mtau, htau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v (mV)
	celsius = 22 (degC)
	dt (ms)
	ena = 55 (mV)
	ek = -70 (mV)
	er = -43 (mV)
	el = -70 (mV)
	gnabar = 0.25 (mho/cm2)    <0,1e9>
	vsna = 0 (mV)
}

STATE {
	m h 
}

ASSIGNED {
	gna (mho/cm2) 
    ina (mA/cm2)
    minf hinf 
    mtau (ms) htau (ms) 
}

LOCAL mexp, hexp 

? currents
BREAKPOINT {
	SOLVE states
    gna = gnabar*(m^3)*h
    ina = gna*(v - ena)
}

UNITSOFF

INITIAL {
    trates(v)
    m = minf
    h = hinf
}

? states
PROCEDURE states() {  :Computes state variables m, h, and n
	trates(v)      :             at the current v and dt.
	m = m + mexp*(minf-m)
	h = h + hexp*(hinf-h)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10
LOCAL qten
? rates
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
LOCAL  alpha, beta, sum

	q10 = 2^((celsius - 22)/10) : was 3^
	qten = 6^((celsius - 22)/10) : was 10^
	
:"m" sodium activation system - JSR
        alpha = -0.36*q10*vtrap((v+49),-3)
        beta =  0.4*q10*vtrap((v+58),20)
        sum = alpha + beta
		mtau = 1/sum
        minf = alpha/sum

:"h" sodium inactivation system - JSR
        alpha = 2.4*q10/(1+exp((v+68-vsna)/3)) + 0.8*qten/(1+exp(v+61.3-vsna))
        beta = 3.6*q10/(1+exp(-(v+21-vsna)/10))
        sum = alpha + beta
		htau = 1/sum
        hinf = alpha/sum




: jsr modified sodium channel - defined in terms of alpha and beta this time
:    am = (0.36*q10*(v+49))/(1-exp(-((v+49)/3)))
:    am = -(0.36*q10*vtrap(-(v+49),3))
:    bm = -(0.40*q10*(v+58))/(1-exp((v+58)/20))
:    bm = (0.40*q10*vtrap((v+58),20))
:    ah = ((2.4*q10)/(1+exp((v+68)/3))) + (0.8*qten/(1+exp(v+61.3)))
:    bh = (3.6*q10)/(1+exp(-(v+21)/10))

:    minf = am/(am+bm)
:    hinf = ah/(ah+bh)
    
:	mtau = 1/(am+bm)
:	htau = 1/(ah+bh)

}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
	TABLE minf, mexp, hinf, hexp DEPEND dt, celsius FROM -100 TO 100 WITH 200

    rates(v)    : not consistently executed from here if usetable_hh == 1
        : so don't expect the tau values to be tracking along with
        : the inf values in hoc

	tinc = -dt * q10
	mexp = 1 - exp(tinc/mtau)
	hexp = 1 - exp(tinc/htau)
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(exp(x/y) - 1)
	}
}

UNITSON
