"""
LELME2150 - Thermal cycles
Homework 1 - Basic cycles

Signature of the function for the Brayton cycle

@author: Antoine Laterre
@date: August 26, 2022
"""

#
#===IMPORT PACKAGES============================================================
#

from CoolProp.CoolProp import PropsSI
from scipy.integrate import quad
from scipy.constants import R
import numpy as np

#
#===BRAYTON CYCLE - TO BE IMPLEMENTED==========================================
#

class Brayton_cycle(object):
   
    
    def __init__(self,inputs,parameters):
        """
        Create a Brayton cycle object.
        """
        p_1,T_1,p_3,T_3         = inputs
        self.p_1,self.T_1       = p_1, T_1                  # pressure [Pa] and temperature [K] at the compressor inlet
        self.p_3,self.T_3       = p_3, T_3                  # pressure [Pa] and temperature [K] at the turbine    inlet
        
        self.eta_pi_c           = parameters['eta_pi_c']    # compressor polytropic efficiency     [-]
        self.eta_pi_t           = parameters['eta_pi_t']    # turbine    polytropic efficiency     [-]
        self.eta_mec_c          = parameters['eta_mec_c']   # compressor mechanical efficiency     [-]
        self.eta_mec_t          = parameters['eta_mec_t']   # turbine    mechanical efficiency     [-]
        self.k_cc               = parameters['k_cc']        # pressure loss coeff. in comb. chamb. [-]
        self.fluid              = parameters['fluid']       # fluid species (e.g. ['N2','O2'])
        self.fluid_prop         = parameters['fluid_prop']  # species proportions (e.g. [0.79, 0.21])
        
    '''def Cp (self, T, f) :
     Cp = PropsSI('CPMASS','P',1e+5,'T',T,'air')
     return Cp
    
    def cp_av(self, Ta, Tb, f) :
        cp_av = 1/(Tb-Ta) * quad(self.Cp,Ta,Tb,args = (f))[0]
        return cp_av 
        '''
 
    def evaluate(self):
        M_molaire_air = 0.79*PropsSI('M','N2') + 0.21*PropsSI('M','O2')
        R_star = R/M_molaire_air
        
     

        # state 1 -> state 2 (compression)

        self.cp_12 = 1005
        gamma_12 = self.cp_12/(self.cp_12-R_star)
        self.p_2 = self.p_3/self.k_cc
        self.T_2 = self.T_1 * ((self.p_2/ self.p_1)** ((gamma_12 - 1) / (gamma_12 * self.eta_pi_c)))

        #self.cp_12 = self.cp_av(self.T_2,self.T_3,0)

        self.h_1 = PropsSI('H', 'T', self.T_1, 'P', self.p_1, 'Air')
        self.h_2 = self.h_1 + self.cp_12*(self.T_2 - self.T_1)
        self.s_1 = PropsSI('S', 'T', self.T_1, 'P', self.p_1, 'Air')
        self.s_2 = self.s_1 + self.cp_12*np.log(self.T_2/self.T_1)-R_star*np.log(self.p_2/self.p_1)

        # state 2 -> state 3 (combustion)
        
        #self.T_3 and self.p_3 are inputs
        self.h_3 = PropsSI('H', 'T', self.T_3, 'P', self.p_3, 'Air')
        self.s_3 = PropsSI('S', 'T', self.T_3, 'P', self.p_3, 'Air')
        

        # state 3 -> state 4 (Expansion)
        self.cp_34 = 1005
        gamma_34 = self.cp_34/(self.cp_34-R_star)
        self.p_4 = self.p_1
        self.T_4 = self.T_3 * ((self.p_4/ self.p_3)** ((gamma_34 - 1) / (gamma_34 * self.eta_pi_t)))
        self.h_4 = self.h_3 + self.cp_34*(self.T_4 - self.T_3)
        self.s_4 = self.s_3 + self.cp_34*np.log(self.T_4/self.T_3)-R_star*np.log(self.p_4/self.p_3)
        # Cycle efficiency - do not modify
        self.w_comp = (self.h_2-self.h_1)/self.eta_mec_c
        self.q_comb = (self.h_3-self.h_2)
        self.w_turb = (self.h_3-self.h_4)*self.eta_mec_t
        self.eta_en = (self.w_turb-self.w_comp)/self.q_comb        
        # Final outputs - do not modify
        self.p = (self.p_1, self.p_2, self.p_3, self.p_4)
        self.T = (self.T_1, self.T_2, self.T_3, self.T_4)
        self.s = (self.s_1, self.s_2, self.s_3, self.s_4)
        self.h = (self.h_1, self.h_2, self.h_3, self.h_4)
        self.states = (self.p,self.T,self.s,self.h)
        
