"""
LELME2150 - Thermal cycles
Homework 2 - Gas turbine

Signature of the function for the gas turbine

@author: Antoine Laterre
@date: September 30, 2022
"""

#
#===IMPORT PACKAGES============================================================
#
from CoolProp.CoolProp import PropsSI
from CoolProp.CoolProp import set_reference_state
from scipy.integrate import quad
from scipy.optimize import fsolve
from numpy import trapz


#
#===BRAYTON CYCLE - TO BE IMPLEMENTED==========================================
#

class gas_turbine(object):
    """
    Class for the simulation of gas turbines

    Parameters
    ----------
    'inputs':
        -> Tuple with the model inputs (SI units = K, Pa).
    'parameters':
        -> Dictionary with the model parameters.
    'display':
        -> Boolean to plot the graphs.
        
    Methods
    -------
    'evaluate':
        -> Compute the cycle
    """
    
    def __init__(self,inputs,parameters,display):
        """
        Create a gas turbine object.
        """
        p_1,T_1,P_e             = inputs
        self.p_1                = p_1                       # [Pa] compressor supply pressure
        self.T_1                = T_1                       # [K]  compressor supply temperature
        self.P_e                = P_e                       # [W]  gas turbine power output
        self.T_3                = parameters['T_3']         # [K]  turbine inlet temperature 
        self.r_c                = parameters['r_c']         # [-]  compression ratio
        self.eta_pi_c           = parameters['eta_pi_c']    # [-]  polytropic efficiency of compressor
        self.eta_pi_t           = parameters['eta_pi_t']    # [-]  polytropic efficiency of turbine
        self.k_cc               = parameters['k_cc']        # [-]  pressure   losses coefficient in combustion chamber
        self.k_mec              = parameters['k_mec']       # [-]  mechanical losses coefficient
        self.air                = parameters['air']         # air  species
        self.air_prop           = parameters['air_prop']    # mol fraction of air  species
        self.alkane             = parameters['alkane']      # alkane composition ['c','h']
        self.display            = display                   # if True, make the plots
    

    T0=0 # [K] reference temperature
    P0=0 # [Pa] reference pressure
    # def cp_species(self, T,comp):
    #         def cp(T) :
    #             return PropsSI('C','T',T,'P',101325,comp)
    #         return cp(T)
    
    def S_species(self, T, P, comp):
            sum=0
            for i in range(len(comp[0])):
                sum+=comp[0][i]*PropsSI('S','T',T,'P',P,comp[1][i])

            return sum
    def H_species(self, T, P, comp):
            sum=0
            for i in range(len(comp[0])):
                sum+=comp[0][i]*PropsSI('H','T',T,'P',P,comp[1][i])
            return sum
    
    def E_species(self, T, P, comp):
        # exergy = H-H0 - T0*(S-S0)
        sum=0
        for i in range(len(comp[0])):
            sum+=comp[0][i]*(PropsSI('H','T',T,'P',P,comp[1][i])-PropsSI('H','T',T0,'P',P0,comp[1][i])-T0*(PropsSI('S','T',T,'P',P,comp[1][i])-PropsSI('S','T',T0,'P',P0,comp[1][i])))
        return sum
    # def cp(self, T1, T2, comp):
    #     #return a function that calculate Cp by integrate Cp(T) from T1 to T2
    #     def func(T2):
    #         result = (quad(lambda x : self.cp_species(x,comp), T1, T2))/(T2-T1)
    #         return result
    #     return func
    
    def R_species(self, comp):
        Mm = comp[0]*28 + comp[1]*32 + comp[2]*44 + comp[3]*18 
        return 1000*8.314/Mm
    
    def cp_species_init(self, T, p, comp):
        result =0
        for i in range(len(comp[0])):
                result+=comp[0][i]*PropsSI('C','T',T,'P',p,comp[1][i])
        return result

    def cp_species(self, T1, T2, p, comp):
        result=0
        for i in range(len(comp[0])):
            T_sat = PropsSI('T','P',p,'Q',0,comp[1][i])
            integrale = quad(lambda x : PropsSI('C','T',x,'P',p,comp[1][i]), T1, min(T2, T_sat))[0] #integrate from T1 to min(T2, T_sat)
            integrale += quad(lambda x : PropsSI('C','T',x,'P',p,comp[1][i]), min(T2, T_sat), T2)[0] #integrate from min(T2, T_sat) to T2
            #print(integrale)
            if(T1<T_sat and T2>T_sat): #if we cross the saturation temperature we need to add the latent heat
                h_liq = PropsSI('H', 'T', T_sat, 'Q', 0, comp[1][i])  
                h_vap = PropsSI('H', 'T', T_sat, 'Q', 1, comp[1][i])
                hlv = h_vap - h_liq
                integrale += hlv
            
            result += comp[0][i]*integrale/(T2-T1)
        return result
    
    def T_compressor_outlet(self, T1, r_c, eta_pi_c, comp):
        T2=T1*(r_c)**(eta_pi_c*self.R_species(comp)/self.cp_species_init(T1, self.p_2, comp))
        cp = self.cp_species(T1, T2, self.p_2, comp)
        return T2
        
    



    def evaluate(self):
        """
        This is the main method of the gas turbine class.
        It evaluates the different thermodynamic quantities at each state of 
        the cycle, as well as some KPI's.
        """
    
        # STATE 1 - INLET OF COMPRESSOR, everything is given in the inputs
        composition = [[self.air_prop[0], self.air_prop[1],0,0],[self.air[0], self.air[1],'CO2','Water']]
        global T0
        global P0
        T0 = self.T_1 # [K] reference temperature
        P0 = self.p_1 # [Pa] reference pressure

        self.s_1 = self.S_species(self.T_1,self.p_1,composition)
        self.h_1 = self.H_species(self.T_1,self.p_1,composition)
        self.e_1 = self.E_species(self.T_1,self.p_1,composition)
        print(self.p_1/1000, self.T_1-273.15, self.h_1/1000, self.s_1/1000, self.e_1/1000)

        # STATE 2 - OUTLET OF COMPRESSOR
        self.p_2 = self.r_c*self.p_1
        R = self.R_species(composition[0])
        #intial guess for T2 and Cp
        cp = self.cp_species_init(self.T_1, self.p_2, composition)
        self.T_2 = self.T_1*(self.r_c)**(self.eta_pi_c*R/cp)

        #iteration to find actual value of cp and T2
        # T2= T1*(r_c)**(eta_pi_c*R/cp)
        # cp = cp_species(T1,T2,p2,composition)
        # ==> fsolve :-) 

        
               
        #print(cp2)



        
        

        # States --------------------------------------------------------------
        self.p           = self.p_1, self.p_2, self.p_3, self.p_4
        self.T           = self.T_1, self.T_2, self.T_3, self.T_4
        self.s           = self.s_1, self.s_2, self.s_3, self.s_4
        self.h           = self.h_1, self.h_2, self.h_3, self.h_4
        self.e           = self.e_1, self.e_2, self.e_3, self.e_4
        self.DAT         = self.p,self.T,self.s,self.h,self.e
        # Combustion paramters ------------------------------------------------
        self.COMBUSTION  = self.LHV,self.e_f,self.excess_air,self.gas,self.gas_prop
        # Mass flow rates -----------------------------------------------------
        self.MASSFLOW    = self.dotm_a,self.dotm_f,self.dotm_g
        # Efficiencies --------------------------------------------------------
        self.ETA         = self.eta_cyclen,self.eta_toten,self.eta_mec,self.eta_cyclex,self.eta_totex,self.eta_rotex,self.eta_combex
        # Energy losses -------------------------------------------------------
        self.DATEN       = self.loss_mec,self.loss_echen
        # Exergy losses -------------------------------------------------------
        self.DATEX       = self.loss_mec,self.loss_rotex,self.loss_combex,self.loss_echex
        # Energy and Exergy pie charts
        if self.display: self.FIG = self.fig_pie_en,self.fig_pie_ex, self.fig_Ts, self.fig_ph
