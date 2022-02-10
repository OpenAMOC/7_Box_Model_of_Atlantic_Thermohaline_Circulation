##Wallis, 2022
#Please visit openamoc.info or investigate the README file to learn more

#import required python modules
import math
from typing import Mapping
from itertools import count
import numpy as np
import matplotlib.pyplot as plt
from openpyxl import Workbook
import sys

#import 2 further python scripts from this package
#Bmodel defines box model functions including rk4 integration, source and sink
import AMOC_Bmodel
#Constants defines constant values such as box boundaries and heat capacity
from AMOC_Constants import *

#create workbook (Excel document with 2 sheets)
print("Using default settings this model takes approx 5 and half minutes to run. Please remain patient.")
wb = Workbook()
ws1 = wb.active
ws1.title = ("Temp, PSU")
wb.save("AMOC_Workbook.xlsx")
ws2 = wb.create_sheet("Time, Flow, Starting values")
print("Starting model...")

#latent heat transport equations
def lht(tc, delta_y):
	'Latent heat transport'
	sat = 6.112*math.exp(17.67*tc/(tc + 243.5))
	dqsdt = (243.5*17.67*0.622*1.0e-3)*sat/(tc + 243.5)**2
	return (1.5*5.1e17*0.8/delta_y)*dqsdt

class Amoc(AMOC_Bmodel.BalanceModel):
	class AmocState(AMOC_Bmodel.State):
		#create variables
		#tos,sos= Temp, Salinity of North Atlantic box
		#tom,som= Temp, Salinity of surface current
		#ton,son= Temp, Salinity of South Atlantic box
		#tod,sod= Temp, Salinity of deep current
		#tas, tam, tan= Temp for atmosphere over north Box, surface current and south box
		_variables_ = ['tos', 'tom', 'ton', 'tod', 'sos', 'som', 'son', 'sod', 'tas', 'tam', 'tan']
		ws1.append(_variables_)

		def __getitem__(self, name):
			if name in ('lht', 'amoc'):
				return getattr(self, name)
			else:
				return super().__getitem__(name)

		@property
		def lht(self): #atmospheric latent heat transfer
			if '_lht' not in self.__dict__:
				tc0 = (self['tam']*(latia[1] - lata[0])+ self['tas']*(lata[1]-latia[1]))/(lata[1] - lata[0])
				lht0 = lht(tc0, ydis[0])
				tc1 = (self['tam']*(latia[2] - lata[2])+ self['tan']*(lata[1]-latia[2]))/(lata[1] - lata[2])
				lht1 = lht(tc1, ydis[1])
				self._lht = np.array([lht0, lht1])
				self.setflags(write=False)
			return self._lht

		@property
		def amoc(self): #flow rate between north and south ocean box
			if '_amoc' not in self.__dict__:
				phi = 1.5264e10*(8.0e-4*(self['son'] - self['sos'])- 1.5e-4*(self['ton'] - self['tos']))
				phi = phi*(phi > 0)
				self._amoc = phi
				self.setflags(write=False)
			return self._amoc

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		
		#set the factors which resist change (inertias)
		#(oceanic heat capacity, ocean box volume, atmospheric heat capacity)
		inertias = [*co, *vo, *ca]
		self.set_inertia(**dict(zip(self.AmocState._variables_, inertias)))
		do_atmos = input("\nInclude atmospheric processes? (The model will function abnormally otherwise) Y/N \n")
		if do_atmos == ("yes") or do_atmos == ("Y") or do_atmos == ("y"):
			do_rad = True #radiative transfer on
			do_srf_hlf = True #heat exchange between the ocean and the atmosphere on
			do_atm_ht = True #atmospheric heat transport on
			do_salt_raise = True #evaporation increases the salinity of the ocean on
			do_salt_amoc = True #amoc salt transport on
			do_heat_amoc = True #amoc heat transport on
			print("\nRunning model...\n")
		else:
			do_rad = False #radiative transfer off
			do_srf_hlf = False #heat exchange between the ocean and the atmosphere  off
			do_atm_ht = False #atmospheric heat transport off
			do_salt_raise = False #evaporation increases the salinity of the ocean off
			do_salt_amoc = True #amoc salt transport on
			do_heat_amoc = True #amoc heat transport on
			print("Running model...")

        #radiative transfer
		if do_rad:
			@self.add_flux(sink='tas')
			def _rad_s(state):
				return areaa[0]*(320.0*(1 - 0.4) - 213.35 - 2.22*state['tas'])

			@self.add_flux(sink='tam')
			def _rad_m(state):
				return areaa[1]*(390.0*(1 - 0.25) - 213.35 - 2.22*state['tam'])

			@self.add_flux(sink='tan')
			def _rad_n(state):
				return areaa[2]*(270.0*(1 - 0.42) - 213.35 - 2.22*state['tan'])

		#heat exchange between ocean and atmosphere
		if do_srf_hlf:
			@self.add_flux(source='tas', sink='tos')
			def _srf_hflx_s(state):
				return areao[0]*(10 - 50*(state['tos'] - state['tas']))

			@self.add_flux(source='tam', sink='tom')
			def _srf_hflx_m(state):
				return areao[1]*(70 - 50*(state['tom'] - state['tam']))

			@self.add_flux(source='tan', sink='ton')
			def _srf_hflx_n(state):
				return areao[2]*(20 - 50*(state['ton'] - state['tan']))

		#atmospheric heat transport
		if do_atm_ht:
			@self.add_flux(source='tam', sink='tas')
			def _aht_m2s(state):
				return (perim[0])*((2.5e13/ydis[0])*(state['tam'] - state['tas'])+ state.lht[0])

			@self.add_flux(source='tam', sink='tan')
			def _aht_m2n(state):
				return (perim[1])*((2.5e13/ydis[1])*(state['tam'] - state['tan'])+ state.lht[1])
	
		#evaporation
		if do_salt_raise:
			@self.add_flux(source='sos', sink='som')
			def _ers_s2m(state):
				return 34.9*perim[0]*(80.0/(360.0*2.5e9))*state.lht[0]

			@self.add_flux(source='son', sink='som')
			def _ers_n2m(state):
				return 34.9*perim[1]*(2.5*80.0/(360.0*2.5e9))*state.lht[1]

		#heat conveyor belt between currents and boxes
		if do_heat_amoc:
			@self.add_flux(sink='tos', source='tod')
			def _heat_amoc_d2s(state):
				return cswt*state.amoc*state['tod']
			
			@self.add_flux(sink='tom', source='tos')
			def _heat_amoc_s2m(state):
				return cswt*state.amoc*state['tos']
            
			@self.add_flux(sink='ton', source='tom')
			def _heat_amoc_m2n(state):
				return cswt*state.amoc*state['tom']

			@self.add_flux(sink='tod', source='ton')
			def _heat_amoc_n2d(state):
				return cswt*state.amoc*state['ton']
		
		#salinity conveyor belt between currents and boxes
		if do_salt_amoc:
			@self.add_flux(sink='sos', source='sod')
			def _salt_amoc_d2s(state):
				return state.amoc*state['sod']

			@self.add_flux(sink='som', source='sos')
			def _salt_amoc_s2m(state):
				return state.amoc*state['sos']

			@self.add_flux(sink='son', source='som')
			def _salt_amoc_m2n(state):
				return state.amoc*state['som']

			@self.add_flux(sink='sod', source='son')
			def _salt_amoc_n2d(state):
				return state.amoc*state['son']

	#printing to results window
	def run(self, ic:Mapping[str, float], nstep, nhist=100, nprint=1000):
		
		time = np.arange(0, nstep, nhist)*self.dt
		
		ic = self.AmocState.from_dict(ic)
		states = self.AmocState(shape=(len(self.AmocState._variables_), nstep//nhist),dtype=float)
		
		states[:, 0] = ic[:]
		for i, state in zip(range(1, nstep), self.iter_states(ic)):
			qut, rmd = divmod(i, nhist)
			if rmd == 0:
				states[:, qut] = state
			if i%nprint == 0:
				print('NSTEP', i)
				print(state)
				#also appends data to workbook each cycle
				test_state = state.tolist()
				ws1.append(test_state)
		return time, states

def main():
	#starting values for tos, tom, ton, tod, sos, som, son, sod, tas, tam, tan
	#give option to use preset starting values
	default_ic = str(input("\nDo you want to use the default starting values for temperature and salinity? Y/N \n"))
	if default_ic == ("yes") or default_ic == ("Y") or default_ic == ("y"):
		ic_values = [4.8, 24.4, 2.7, 2.7, 34.4, 35.6, 34.9, 34.9, 4.7, 3.3, 0.9]
		print("\nDefault values are: "+ str(ic_values) +"\n")
	else:
		print("\nPlease enter starting values\n")
		tos_q = float(input("Starting value for tos: "))
		tom_q = float(input("Starting value for tom: "))
		ton_q = float(input("Starting value for ton: "))
		tod_q = float(input("Starting value for tod: "))
		sos_q = float(input("Starting value for sos: "))
		som_q = float(input("Starting value for som: "))
		son_q = float(input("Starting value for son: "))
		sod_q = float(input("Starting value for sod: "))
		tas_q = float(input("Starting value for tas: "))
		tam_q = float(input("Starting value for tam: "))
		tan_q = float(input("Starting value for tan: "))
		ic_values = [tos_q, tom_q, ton_q, tod_q, sos_q, som_q, son_q, sod_q, tas_q, tam_q, tan_q,]
		print("\nYour chosen values are: " +str(ic_values) +"\n")

	ic = dict(zip(Amoc.AmocState._variables_, ic_values))
	NSEC_YR = 86400*365
	dt_yr = 0.01
	
	#integration method of core differential equation
	dt, nstep = dt_yr*NSEC_YR, int(5.0e5)
	scheme = 'rk4'
	
	#run until the equilibrium has been reached
	model = Amoc(scheme=scheme, dt=dt)
	time0, states0 = model.run(ic=ic, nstep=nstep)
	
	#freshwater hosing experiment in northern box
	print('\n\nFirst run complete\n\n')
	print('\n\nInitiating salinity drop experiment\n\n')
	ws1.append(["First run complete, below is the salinity drop experiment"])
	ic1 = np.array(states0[:, -1])
	ic1[Amoc.AmocState._variables_.index('son')] -= 0.7
	ic1 = dict(zip(Amoc.AmocState._variables_, ic1))
	time1, states1 = model.run(ic=ic1, nstep=nstep)
	
	time1 += time0[-1]
	
	for time in (time0, time1):
		time /= NSEC_YR # in year
		
	#for legend
	def txt2tex(txt):
		return f'${txt[0].capitalize()}_{txt[2].capitalize()}^{txt[1]}$'
	
	#colours for different regions
	colors = [plt.get_cmap("tab10")(i) for i in range(4)]
	
	#the first 3000 years
	def xy0(v):
		nstep = 3000
		return time0[:nstep], states0[v][:nstep]
	
	#the whole 10000 years
	def xy1(v):
		return map(np.concatenate,((time0, time1), (states0[v], states1[v])))
		
	#plot variables across 3 axes
	for j, xy in enumerate((xy0, xy1)):
		fig, (ax0, ax1, ax2) = plt.subplots(ncols=1, nrows=3, figsize=(8, 8))
		for i, v in enumerate(('tos', 'tom', 'ton', 'tod')):
			ax0.plot(*xy(v), label=txt2tex(v),linestyle='--', color=colors[i])
		for i, v in enumerate(('tas', 'tam', 'tan')):
			ax0.plot(*xy(v), label=txt2tex(v), color=colors[i])
		for i, v in enumerate(('sos', 'som', 'son', 'sod')):
			ax1.plot(*xy(v), label=txt2tex(v), color=colors[i])
			
		time, amoc = xy('amoc')
		ax2.plot(time, amoc*1.0e-6, color=colors[0])
		ax0.set_ylabel(r'Temperature (C$^\circ$)')
		ax1.set_ylabel(r'Salinity (psu)')
		ax2.set_ylabel(r'$\Phi$ (10$^6$ Sv)')
		ax2.set_xlabel(r'Time (yr)')
		
		for ax in (ax0, ax1):
			ax.legend(ncol=2)
			ax.grid()
		
		ax2.grid()
		#add data to workbook
		time_list = time.tolist()
		ws2.append(time_list)
		flowrate = (amoc*1.0e-6)
		flowrate_list = flowrate.tolist()
		ws2.append(flowrate_list)
		ws2.append(ic_values)
		#save workbook
		print("\nDo Not Close Program Until Prompted. Please close graph window when ready.\n")
		wb.save("AMOC_Workbook.xlsx")
		#save and show figures
		plt.tight_layout()
		plt.savefig(f'timeseries{j}.png')
		plt.show()
		plt.close(fig)

if __name__ == '__main__':
    main()

#sign off program
print("Items saved, you may now close the program.")
