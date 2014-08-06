""" Weighted Histogram Analysis Method (WHAM). Optimized way to analyse
    free energy and free energy curves from umbrella sampling simulations.
    The input are the coupling parameters for each intermediate states
    and the associated initial-final state potential energy gaps.

    Options to change the convergence criteria, the number of datapoints
    used for convergence are included. 

    Functions to plot and analyze the gaps as well as the weights 
    for each sampling window are available. 

    The maxima and minima of the initial-final state potential energy
    gaps defines the reaction coordinate used for plotting, with a
    default resolution of 0.1 eV.

"""

from ase.units import kB
import numpy as np
import pylab as pl
from printer import *

class WHAM:

    def __init__(self, energy_gaps, coupling_para, T = 300, res = 0.1, 
                 iter = 100, npt = None, free_w = None, conv=1e-4, pr=None):

        self.energy_gaps = energy_gaps # Energy gaps. array((w,Nw))
                             # w is the number of sampling windows
                             # Nw is the number of datapoints per window
                             # Make sure Nw is the same for each w
        self.coupling_para = coupling_para  # The umbrella sampling coupling
                      # parameters for each window w. array(0,eta1,eta2...1)
        if isinstance(T, float):
	    self.T = np.array([T])
	else:
            self.T = T        # Simulation temperature, default 300
        self.res = res    # data is condensed to arange(dE_min,dE_max,res)
        self.iter = iter  # number of iterations for the free energies 
                          # of each window, fj.
        self.npt = npt    # number between each datapoint, so instead of  
                          # using say 1.000.000 points one can quickly run
                          # 1.000.000 / npt points (good for initial fj)
        self.free_w = free_w # free energies of perturbation for each 
                             # sampling window. If none an initial guess
                             # is provided. array(0,0.1,0.2,0.3,0.4)
        self.no = None    # number of point per window
        self.conv = conv  # Convergence criteria

        self.cycle = False # have the free_w been solved for?

        self.Pu = None      # hold on to unbiased probs when calc
        self.Pt = None      # hold on to total probability
        self.Ft = None      # hold on to free energy curve
        self.Pb = None      # hold on to biased probability curves
        self.pl_arr = None  # hold on to reaction coordinate array

        self.weights = None # hold on to weights when calculated

        self.plot = None

        self.pr=free_printer(pr)
        self.pr.header()
        self.update_params()


    def get_free_energies(self, it=None, pr=None, conv=None):
        """ solve for the free energies of perturbation for each 
            window.
        """
        self.pr.sub_header('get_free_energies()')
        dE = self.energy_gaps
        cp = self.coupling_para
        free = self.free_w
        no = self.no
        w = len(cp)
        nt = self.npt

        if it is None:
            it = self.iter
 
        if conv is None:
            conv = self.conv # Default value

        # Already made sure that dE and eta are in order
        # Initialize
	if len(self.T)==1:
		betaval = 1./(self.T * kB)
		beta = np.zeros_like(cp)+betaval
	else:
		beta = 1./(self.T * kB)
        ebf2 = np.zeros(w)
        ebf = np.zeros(w)
        fact = np.zeros(w)

        ebf[:] = np.exp(beta*free)
        fact[:] = no * np.exp(beta*free)

        # Start the iteration
        for I in range(it):

            old = np.zeros(w)
            new = np.zeros(w)

            if pr != None:
                self.pr.f_print('Iteration no. : '+str(I)+' Free energies : '+str(1./beta * np.log(ebf)))
                #print 'Iteration no. :', I, 'Free energies :', \
                 #                           1./beta * np.log(ebf)

            old[:] = ebf

            for k in range(w):
                ebfk = 0

                for i in range(w):

                    bottom = (np.exp(-beta*cp*dE[i,::nt].reshape(len(dE[i,::nt]),-1))\
                              *fact).sum(axis=1)
                    ebfk += (np.exp(-beta[k]*cp[k]*dE[i,::nt]) / bottom).sum()

                ebf2[k] = ebfk
                ebf[k]  = 1./(ebf[0]*ebfk)
                fact[k] = no * ebf[k]
                new[k] += ebf[k]

            if abs((1./beta * np.log(old) - 1./beta * np.log(new)).sum()) \
                    <= conv:
                self.pr.f_print('converged : '+str(abs((1./beta * np.log(old) - \
                                          1./beta * np.log(new)).sum()) <= conv)+ \
                                     ' Norm. free energies : '+ \
                                     str(1./beta * np.log(new) - 1./beta * np.log(new[0])))
                #print 'converged :', abs((1./beta * np.log(old) - \
                #                          1./beta * np.log(new)).sum()) <= conv, \
                #                     'Norm. free energies :', \
                #                     1./beta * np.log(new) - 1./beta * np.log(new[0])
                break    

        free_new = 1./beta * np.log(ebf)

        # Normalize
        free_new[:] -= free_new[0]
        self.free_w = free_new

        # Have the free energies been determined?
        self.cycle = True 
        self.pr.sub_footer('get_free_energies()')

    def get_biased_probabilities(self, res=None, free=False):
        # Run over all energy gaps in the arange(min,max,res)
        # and collect probabilities from delta function
        # If free is not none the free energies of these probabilites
        # is returned as well.
        self.pr.sub_header('get_biased_probabilities()')
        dE = self.energy_gaps
        cp = self.coupling_para
        nt = self.npt
        no = self.no
        win = len(cp)

        beta = 1./(self.T * kB)

        if res is None:
            res = self.res

        # Make sampling array
        x = np.arange(dE.min(),dE.max(),res)

        # Collect prob
        P = np.zeros((win,len(x)))

        for w in range(win):
            for a in range(len(x)):
                x1 = dE[w, ::nt] >= x[a]
                x2 = dE[w, ::nt] <  x[a] + res
                x3 = x1 == x2
                P[w,a] += float(x3.sum()) / no

        self.pl_arr = x
        self.Pb = P

        if free:
            return x, -1./beta * np.log(P)
        self.pr.sub_footer('get_biased_probabilities()')

    def get_unbiased_probabilities(self, res=None, free=False):
        # Assume free_w have been solved for? or say:
        # if self.cycle:
        #     fj = self.free
        # else:
        #     run free energy thingy in silent mode

        self.pr.sub_header('get_unbiased_probabilities()')
        dE = self.energy_gaps
        cp = self.coupling_para
        fj = self.free_w
        nt = self.npt
        no = self.no
        win = len(cp)

	if len(self.T)==1:
		betaval = 1./(self.T * kB)
		beta = np.zeros_like(cp)+betaval
	else:
		beta = 1./(self.T * kB)

        if res is None:
            res = self.res

        # Make sampling array
        x = np.arange(dE.min(),dE.max(),res)

        # Collect prob
        P = np.zeros((win,len(x)))

        for w in range(win):
            for a in range(len(x)):
                x1 = dE[w, ::nt] >= x[a]
                x2 = dE[w, ::nt] <  x[a] + res
                x3 = x1 == x2
                
                # Apply weights
                x4 = x3 * np.exp(beta[w]*(cp[w]*dE[w,::nt] - fj[w]))

                # sum up
                P[w,a] += x4.sum() / no

        self.Pu = P
        self.pl_arr = x

        if free:
            return x, -1./beta * np.log(P)

        self.pr.sub_footer('get_unbiased_probabilities()')
    def get_total_probability_curve(self, res=None, free=False):
        # Combine weights and unbiased probabilities
        # save the object to total prob curve
        self.pr.sub_header('get_total_probability_curve()')

        if res is None:
           res = self.res

        beta = 1./(self.T * kB)        

        # Assume free energies are solved for?

        # This needs a cleaner write up... need to make sure that
        # the p_x and Pu arrays have the same length...
        # Maybe remove the resolution option for all but the intitial
        # call, and only have the option via update_params...
        if self.weights is None:
            self.get_weights(res=res)

        if self.Pu is None:
            self.get_unbiased_probabilities(res=res)
        elif  len(self.Pu[0,:]) != len(self.weights[0,:]):
            self.get_unbaised_probabilities(res=res)
            
        self.Pt = (self.weights * self.Pu).sum(axis=0)
        self.Ft = -1./beta[0] * np.log(self.Pt)

        if free:
            return self.pl_arr, self.Ft
        self.pr.sub_footer('get_total_probability_curve()')

    def get_weights(self, res=None):
        # Collect the weights for each sampling window

        # Assume free energies are solved for? or 
        # run free energy thingy in silent mode
        self.pr.sub_header('get_weights()')

        dE = self.energy_gaps
        cp = self.coupling_para
        fj = self.free_w
        nt = self.npt
        no = self.no
        win = len(cp)
        
        if res is None:
            res = self.res

	if len(self.T)==1:
		betaval = 1./(self.T * kB)
		beta = np.zeros_like(cp)+betaval
	else:
		beta = 1./(self.T * kB)

        # Make sampling array
        x = np.arange(dE.min(),dE.max(),res)

        # Collect weights
        p_x = np.zeros((win,len(x)))

        for w in range(win):
            i = 0
            for a in x:
                p_x[w, i] += no*np.exp(-beta[w] * (cp[w] * a - fj[w])) \
                             /(no*np.exp(-beta * (cp * a - fj))).sum()
                i += 1
        
        self.pl_arr = x
        self.weights = p_x
        self.pr.sub_footer('get_weights()')

    def get_free_energy_curve(self, res=None):
        self.pr.sub_header('get_free_energy_curve()')
       
        x = self.pl_arr
        y = self.Ft 
	y2 = y+x
        
        ind = len(x)/10

        c1 = np.polyfit(x[ind:-1*ind], y[ind:-1*ind], 2)
        c2 = np.polyfit(x[ind:-1*ind], y2[ind:-1*ind], 2)

	self.pr.f_print('Fitting parameters for curve 1: '+str(c1))
	self.pr.f_print('Fitting parameters for curve 2: '+str(c2))
        
	#print c1
	#print c2
        #print ind

        yfit = c1[0]*x**2 + c1[1]*x + c1[2]
        y2fit = c2[0]*x**2 + c2[1]*x + c2[2]

        dU = np.min(y2)-np.min(y)
        dUfit = np.min(y2fit)-np.min(yfit)
        #print np.min(y)
        #print np.min(y2)
        #print np.min(yfit)
        #print np.min(y2fit)
        xmin = np.argmin(y)
        xminfit = np.argmin(yfit)
        x2min = np.argmin(y2)
        x2minfit = np.argmin(y2fit)

        lamb = y[x2min] - y[xmin] 
        lambfit = yfit[x2minfit] - yfit[xminfit] 
        lamb2 = y2[xmin] - y2[x2min] 
        lamb2fit = y2fit[xminfit] - y2fit[x2minfit] 

        self.pr.f_print('dU determined from data points: '+str(dU))
        self.pr.f_print('dU determined from fitted curves: '+str(dUfit))
        #print dU
        #print dUfit

        self.pr.f_print('lambda determined from data points (side 1): '+str(lamb))
        self.pr.f_print('lambda determined from fitted curves (side 1): '+str(lambfit))
        self.pr.f_print('lambda determined from data points (side 2): '+str(lamb2))
        self.pr.f_print('lambda determined from fitted curves (side 2): '+str(lamb2fit))
        #print lamb
        #print lambfit
        #print lamb2
        #print lamb2fit

        bar = (dU - lamb)**2 / (4*lamb)
        barfit = (dUfit - lambfit)**2 / (4*lambfit)

        self.pr.f_print('barrier height determined from data points (side 1): '+str(bar))
        self.pr.f_print('barrier height determined from fitted curves (side 1): '+str(barfit))
        #print bar
        #print barfit
	pl.plot(x, y, color='red', marker='x')
	pl.plot(x, yfit, color='red', linestyle='-', linewidth=3)
	pl.plot(x, y2, color='blue', marker='x')
	pl.plot(x, y2fit, color='blue', linestyle='-', linewidth=3)
	pl.ylabel(r'Free Energy [$e$V]', fontsize=24)
	pl.xlabel(r'Reaction Coordinate, $\Delta E$ [$e$V]',fontsize=24)
	pl.xticks(fontsize=24)
	pl.yticks(fontsize=24)
	self.plot = pl

	#tempfil = open('tempfil.txt', 'w')
	#for n in range(len(x)):
	#	tempfil.write(str(x[n])+'\t'+str(y[n])+'\t'+str(y2[n])+'\n')
	#tempfil.close()

        self.pr.sub_footer('get_free_energy_curve()')
    def get_gap_distribution(self, res=None):
        print 'distr'

    def update_params(self, npt=None):
	if npt:
	    self.npt=npt
        # Make sure npt is in order
        if self.npt is None:
            self.npt = 1
        # Make sure it is an integer
        else: 
            self.npt = int(self.npt)

        # Make initial free_w guess if none is provided
        if self.free_w is None:
            free = np.zeros(len(self.coupling_para))
            U = (self.energy_gaps[-1] + self.energy_gaps[0]).mean()/2.0
            lamb = abs((self.energy_gaps[-1] - self.energy_gaps[0]).mean()/2.0)
            free[0]=0
            free[-1]=U
            if abs(U) <= 0.1:
                meanval = (len(self.coupling_para)+1)/2.0
                for n in range(1,len(self.coupling_para)-1):
                    if n < meanval:
                        free[n] = (lamb/4.0)*self.coupling_para[n]
                    elif n > meanval:
                        free[n] = (lamb/4.0)*(1-self.coupling_para[n])+U
                    elif n == meanval:
                        free[n] = (lamb/4.0)
                        
            self.free_w = free

        # determine the number of dpts per window
        self.no = len(self.energy_gaps[0,:]) / self.npt
        #print self.no

    def show_plot(self, save=None):
        if self.plot:
	    if not save:
                self.plot.show()
            else:
                self.plot.savefig(save)
