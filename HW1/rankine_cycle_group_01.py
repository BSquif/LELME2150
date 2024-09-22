"""
LELME2150 - Thermal cycles
Homework 1 - Basic cycles

Signature of the function for the Rankine cycle

@author: Antoine Laterre
@date: August 26, 2022
"""

#
#===IMPORT PACKAGES============================================================
#

from CoolProp.CoolProp import PropsSI
# you can add other basic scientific packages (SciPy, NumPy, ...)

#
#===RANKINE CYCLE - TO BE IMPLEMENTED==========================================
#

class Rankine_cycle(object):
    """
    Class for the simulation of the Rankine cycle

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (SI units = K, Pa).
    'parameters':
        -> Dictionary with the model parameters.
        
    Methods
    -------
    'evaluate':
        -> Compute the cycle and evaluate some KPI's.
    """
    
    def __init__(self,inputs,parameters):
        """
        Create a Rankine cycle object.
        """
        T_1,p_3,T_3,P_el,LHV = inputs
        self.T_1                = T_1                       # temperature [K] at the pump inlet
        self.p_3,self.T_3       = p_3, T_3                  # pressure [Pa] and temperature [K] at the turbine inlet
        self.P_el, self.LHV     = P_el, LHV                 # demanded power [W] and fuel LHV [J/kg]
        
        self.eta_pump           = parameters['eta_pump']    # pump internal efficiency          [-]
        self.eta_is_t           = parameters['eta_is_t']    # turbine isentropic efficiency     [-]
        self.eta_mec_t          = parameters['eta_mec_t']   # turbine mechanical efficiency     [-]
        self.eta_gen            = parameters['eta_gen']     # steam generator efficiency        [-]
        self.k_gen              = parameters['k_gen']       # pressure losses coeff. in s. gen. [-]
        self.fluid              = parameters['fluid']       # fluid (e.g. 'H2O')

    # def cpComputation(T1, T2, p1, p2, fluid)
    #     # compute your cp by using a constant pressure that is equal to the mean value of pressure at temperature 1 and pressure at temperature 2
    #     cp = 
        
    def evaluate(self):
        """
        This is the main method of the Rankine cycle class.
        It evaluates the different thermodynamic quantities at each state of 
        the Rankine cycle, as well as some KPI's.
        """
    
        # >>>>>             <<<<< #    
        # Replace with your model #
        # >>>>>             <<<<< # 
        
        #inital data :
        #T_1, p_3, T_3, P_el, LHV
        #fluid = 'H2O'
        #eta_pump, eta_is_t, eta_mec_t, eta_gen, k_gen

        #state 1
        self.p_1 = PropsSI('P','T', self.T_1 ,'Q',0,self.fluid) #+3K because T_1 is subcooled
        self.h_1 = PropsSI('H','T', self.T_1,'Q',0,self.fluid)
        self.s_1 = PropsSI('S','T', self.T_1,'Q',0,self.fluid)
        #self.x_1 = PropsSI('Q','T', self.T_1,'P',self.p_1,self.fluid) #x_1 = 0
        self.x_1 = 0
        #print(PropsSI('Q','T', self.T_1,'P',self.p_1,self.fluid))

        #state 1 -> state 2 (compression, pump is represented with the internal efficiency model)
        self.p_2 = self.p_3 / self.k_gen #pressure loss in the steam generator
        delta_T  = ((PropsSI('V','T',self.T_1,'P',self.p_1,self.fluid) * (self.p_2 - self.p_1)) / PropsSI('CPMASS','T',self.T_1,'P',(self.p_1 + self.p_2)/2,self.fluid)) * (1/self.eta_pump -1)
        self.T_2 = self.T_1 + delta_T
   
        self.h_2 = PropsSI('H','T',self.T_2,'P',self.p_2,self.fluid)
        self.s_2 = PropsSI('S','T',self.T_2,'P',self.p_2,self.fluid)
        self.x_2 = PropsSI('Q','T',self.T_2,'P',self.p_2,self.fluid) #x_2 = [-]
        #print(self.x_2)

        #state 2 -> state 3 (steam generator)
        self.h_3 = PropsSI('H','T',self.T_3,'P',self.p_3,self.fluid)
        self.s_3 = PropsSI('S','T',self.T_3,'P',self.p_3,self.fluid)
        self.x_3 = PropsSI('Q','T',self.T_3,'P',self.p_3,self.fluid) #x_3 = [-]
        #print(self.x_3)

        #state 3 -> state 4 (expansion, Adiabatic expansion in the turbine
        #down to the saturation pressure (use the isentropic model),
        #The mechanical efficiency of the turbine is taken into account)
        self.T_4 = self.T_1 + 3
        s4s = self.s_3
        x4s= (s4s - PropsSI('S','T',self.T_4,'Q',0,self.fluid)) / (PropsSI('S','T',self.T_4,'Q',1,self.fluid)-PropsSI('S','T',self.T_4,'Q',0,self.fluid))
        h4s = PropsSI('H','T',self.T_4,'Q',x4s,self.fluid)
        h4 = self.h_3 - self.eta_is_t*(self.h_3-h4s)
        self.p_4 = PropsSI('P','T',self.T_4,'Q',x4s,self.fluid)
        self.h_4 = h4
        self.s_4 = PropsSI('S','T',self.T_4,'Q',x4s,self.fluid)
        self.x_4 = x4s
       
        #print(self.x_4)
        #print(PropsSI('Q','T',300,'P',100000,self.fluid))

        # ...
        # ...
        # ...
        
        # Cycle efficiency - do not modify
        self.W_pump = (self.h_2-self.h_1)
        self.Q_comb = (self.h_3-self.h_2)/self.eta_gen
        self.W_turb = (self.h_3-self.h_4)*self.eta_mec_t
        self.eta_en = (self.W_turb-self.W_pump)/self.Q_comb
        self.dot_m_f= self.P_el/self.LHV/self.eta_en
        # Final outputs - do not modify
        self.p = (self.p_1, self.p_2, self.p_3, self.p_4)
        self.T = (self.T_1, self.T_2, self.T_3, self.T_4)
        self.s = (self.s_1, self.s_2, self.s_3, self.s_4)
        self.h = (self.h_1, self.h_2, self.h_3, self.h_4)
        self.x = (self.x_1, self.x_2, self.x_3, self.x_4)
        self.states = (self.p,self.T,self.s,self.h,self.x)

