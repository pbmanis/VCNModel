# calyx_biophysics.hoc
# The Calyx Biophysics (Channel insertion)
# Pulled from Calyx.hoc just to keep the code cleaner
#
# Every different kind of section has it's own conductance levels
# and can have different channels
#
# Paul B. Manis, Ph.D.
# 25 Sept. 2007
# Modified to use gca for HH formulation of calcium current
# 14 Oct 2007
# build the biophysics... 
def biophys(mode = 0): 
    for sec in stalk:
     	if(mode == 1):
     		sec.insert(klt)
     		sec.klt.ek_klt = Ek
    		sec.insert(kht)
    		sec.kht.ek_kht = Ek
    		sec.insert(na)
    		sec.na.ena_na = 50
    		sec.insert(ih)
    		sec.ih.eh_ih=-43
    		sec.insert(CaPCalyx)
    		sec.CaPCalyx.eca_CaPCalyx=43.9

      	sec.gcapbar_CaPCalyx = gca_st 
     	sec.gnabar_na = gna_st
    	sec.gkhtbar_kht = ghvk_st
    	sec.gkltbar_klt = glvk_st
    	sec.ghbar_ih = gh_st
    
    print 'inserted biophysics'
    return
    
#     forsec branch {
#   if(mode == 1) {
#       insert klt ek_klt = Ek
#       insert kht ek_kht = Ek
#       insert na ena_na = 50
#       insert ih eh_ih=-43
#       insert CaPCalyx eca_CaPCalyx=43.9
#       }
#       gcapbar_CaPCalyx = gca_br 
#   gnabar_na = gna_br
#   gkhtbar_kht = ghvk_br
#   gkltbar_klt = glvk_br
#   ghbar_ih = gh_br
#   }
#    forsec neck {
#   if(mode == 1) {
#       insert klt ek_klt = Ek
#       insert kht ek_kht = Ek
#       insert na ena_na = 50
#       insert ih eh_ih=-43
#       insert CaPCalyx eca_CaPCalyx=43.9
#       }
#       gcapbar_CaPCalyx = gca_nk 
#   gnabar_na = gna_nk
#   gkhtbar_kht = ghvk_nk
#   gkltbar_klt = glvk_nk
#   ghbar_ih = gh_nk
#   }
#    forsec swelling {
#   if(mode == 1) {
#       insert klt ek_klt = Ek
#       insert kht ek_kht = Ek
#       insert na ena_na = 50
#       insert ih eh_ih=-43
#       insert CaPCalyx eca_CaPCalyx=43.9
#       }
#       gcapbar_CaPCalyx = gca_sw 
#   gnabar_na = gna_sw
#   gkhtbar_kht = ghvk_sw
#   gkltbar_klt = glvk_sw
#   ghbar_ih = gh_sw
#   }
#   forsec tip{
#   if(mode == 1) {
#       insert klt ek_klt = Ek
#       insert kht ek_kht = Ek
#       insert na ena_na = 50
#       insert ih eh_ih=-43
#       insert CaPCalyx eca_CaPCalyx=43.9
#   }
#   gcapbar_CaPCalyx = gca_tp
#   gnabar_na = gna_tp
#   gkhtbar_kht = ghvk_tp
#   gkltbar_klt = glvk_tp
#   ghbar_ih = gh_tp
#   }
# 
# }
