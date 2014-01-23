// calyx_electrodes.hoc
//
// insert and manage the vc and cc electrodes. 
// note that the default setup puts the model in current clamp
// use choose_cclamp and choose_vcclamp below to change the recording mode
// Pulled from calyx5.hoc, 11/5/2007
// Paul B. Manis, Ph.D.
// UNC Chapel Hill
//

proc insert_clamps() {
	if(GBCFLAG == 1) {
		printf("Iclamp in GBC Soma\n")
		soma {
		    icgbc = new IClamp(0.5)
			icgbc.amp = icAmpDef
			icgbc.del = 2
			icgbc.dur = icDurDef
		}
	} 
	{
		printf("Iclamp in Cut Calyx Axon[1]\n")
		axon[axonnode] {
		    ic = new IClamp(iclocation)
			ic.amp = icAmpDef
			ic.del = 2
			ic.dur = icDurDef
		}
	}
	saveic()

	axon[axonnode] {
		vc = new SEClamp(vclocation)
		vc.dur1 = 2
		vc.dur2 = 10
		vc.dur3 = 5
		vc.amp1 = -65
		vc.amp2 = -75
		vc.amp3 = -65
		vc.rs = 1e9// make Rs small to turn on clamp 
		savevc()
	}

}

proc setdeficlamp() {
	if(GBCFLAG == 1) {
		soma {
			icgbc.loc = iclocation
			icgbc.amp = icAmpDef
			icgbc.del = 2
			icgbc.dur = icDurDef
		} 
	} 
	{
		axon[axonnode] {
			ic.amp = icAmpDef
			ic.del = 2
			ic.dur = icDurDef
		}
	}
}

proc setdefvclamp() {
	axon[axonnode] {
		vc.dur1 = 2
		vc.dur2 = 10
		vc.dur3 = 5
		vc.amp1 = -80
		vc.amp2 = -85
		vc.amp3 = -80
	}
}

		
proc savevc() {
 if(vsaveflag == 0) {
 	vcdur1 = vc.dur1
	vcdur2 = vc.dur2
	vcdur3 = vc.dur3
	vcamp1 = vc.amp1
	vcamp2 = vc.amp2
	vcamp3 = vc.amp3
	vcrs = vc.rs
	vcsaveflag = 1
 }
 
}
proc restorevc() {
	vc.dur1 = vcdur1
	vc.dur2 = vcdur2
	vc.dur3 = vcdur3
	vc.amp1 = vcamp1
	vc.amp2 = vcamp2
	vc.amp3 = vcamp3
	vc.rs = vc.rs
	vcsaveflag = 0
}

proc saveic() {
	if(icsaveflag == 0) {
		if(GBCFLAG == 1) {
			soma {
				icgbcamp = ic.amp
				icgbcdur = ic.dur
				icgbcdel = ic.del
			}
		}
		{
			axon[axonnode] {
				icamp = ic.amp
				icdur = ic.dur
				icdel = ic.del
			}
		}				
		icsaveflag = 1
	}
}
proc restoreic() {
	if(GBCFLAG == 1) {
		soma {
			icgbc.amp = icgbcamp
			icgbc.dur = icgbcdur
			icgbc.del = icgbcdel
		}
	} 
	 {
		axon[axonnode] {
			ic.amp = icamp
			ic.dur = icdur
			ic.del = icdel
		}
	}
		
	icsaveflag = 0
}


proc choose_iclamp() {
	if(GBCFLAG == 1) {
		soma {	restoreic() }
	} 
	{
		axon[axonnode] { restoreic() }
	}
	clamp_mode = 0
	clampstr = "CC"
	axon[axonnode] {
		vc.rs = 1e9
	}
}

proc choose_vclamp() {
	axon[axonnode] {
		restorevc()
		vc.rs = 0.0001 // and... this too 
		clamp_mode = 1
		clampstr = "VC"
	}
	if(GBCFLAG == 1) {
		soma {
			saveic()
			icgbc.del = 1e9 // and change IC
			icgbc.amp = 0
		}	 
	}
	{
		axon[axonnode] {
			saveic()
			ic.del = 1e9 // and change IC
			ic.amp = 0
		}
	}


}


