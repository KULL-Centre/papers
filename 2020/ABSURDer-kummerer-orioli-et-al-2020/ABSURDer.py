import numpy as np
import scipy.optimize as optimize
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import scipy.stats as scs
from itertools import combinations
import warnings

# rcParams for plots
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['axes.labelsize']  = 18
mpl.rcParams['legend.fontsize'] = 14
warnings.filterwarnings('ignore')

class ABSURDer:

    def __init__( self, rex, rmd, eex = '', out = 'results', thetas = np.array([100,1000,10000]), idx = [], methyl_list = None ):

        """
        Class constructor

        Parameters
        ----------
        rex : str or numpy.ndarray
            path to or numpy ndarray with dimensions nrates*nmethyls containing the experimental rates.
        rmd : str or numpy.ndarray
            path to or numpy ndarray with dimensions nrates*nmethyls*nblocks containing the simulated rates.
        eex : str or numpy.ndarray
            path to or numpy ndarray with dimensions nrates*nmethyls containing the experimental errors.
            If empty, a toy model will be built from rex.
            Default: ''.
        out : str
            path to a pickle file where to store results.
            Default: 'results'.
        thetas : numpy.array
            array of fudge parameters to be used for reweighting.
            Default: np.array([0,100,1000,10000]).
        idx : list
            list of methyls to exclude from the analysis.
            Default: [].
        methyl_list: str
            path to the list of methyl group names. The ordering is expected to be the same as in the rate matrices.
        """

        self.rex = self.load_rates( rex, "experimental" )    # experimental rates
        self.rmd = self.load_rates( rmd, "simulated" )       # simulated rates
        if eex == '' and len(self.rex.shape) == 3:
            self.build_toyexp()                              # no errors provided: build a toy model from rex
            print("# No experimental errors provided: toy model built.")
        elif eex == '' and len(self.rex.shape) != 3:         # no errors and no blocks provided
            raise ValueError('Experimental errors not provided, but the size of rex is not compatible with the construction of a toy model.' )
        else:
            self.eex = self.load_rates( eex, "errors on" )   # experimental errors
        self.idx = idx

        # consistency checks
        if len(self.rex.shape) != 2:
            raise ValueError("rex dimension has to be nrates x nmethyls")
        if len(self.eex.shape) != 2:
            raise ValueError("eex dimension has to be nrates x nmethyls")
        if len(self.rmd.shape) != 3:
            raise ValueError("rmd dimension has to be nrates x nmethyls x nblocks")
        if self.rex.shape[0] != self.eex.shape[0] or self.rex.shape[0] != self.rmd.shape[0]:
            raise ValueError("The number of rates must be identical in rex, eex and rmd")
        if self.rex.shape[1] != self.eex.shape[1] or self.rex.shape[1] != self.rmd.shape[1]:
            raise ValueError("The number of methyl groups must be identical in rex, eex and rmd")

        if methyl_list != None:
            self.load_methyl_list(methyl_list)

        if len(self.idx) > 0:
            self.ignore_methyls()                                    # ignore some methyls in the analysis

        self.ths = -np.sort(-thetas)                                 # set of fudge parameters
        self.out = out                                               # path to file where to save results
        self.r   = self.rex.shape[0]                                 # number of rates
        self.m   = self.rex.shape[1]                                 # number of methyls
        self.b   = self.rmd.shape[-1]                                # number of blocks
        self.emd = np.std( self.rmd, axis = -1 ) / np.sqrt( self.b ) # errors on simulated rates
        self.rav = np.mean( self.rmd, axis = -1 )                    # average simulated rates
        self.w0  = np.full( self.b, 1./self.b )                      # initial weights
        self.ix2 = []                                                # initial chi2 before optimization
        self.res = {}                                                # dictionary of optimized weights for each theta
        self.phi = []                                                # fraction of effective frames
        self.chi = []                                                # reduced chi squared corresponding to optimized weights
        self.cor = []                                                # Pearson correlation corresponding to optimized weights

        self.specdens_load = False                                   # spectral densities have been loaded
        self.rot_load      = False                                   # rotamers have been loaded

        for r in range(self.r):
            self.ix2.append( self.chi2r(r, self.w0) )
        self.ix2.append( self.chi2r(-1, self.w0) )

        print("\n# INFO ON THE DATASET")
        print(f"# Number of methyls:  {self.m}")
        print(f"# Number of rates:    {self.r}")
        print(f"# Number of blocks:   {self.b}")
        print(f"# Overall chi square: {self.ix2[-1]:.2f}" )
    #----------------------------------------------------------------------------------------------------------------

    def load_rates( self, r, tp = "" ):

        """
        Loads the provided rates depending on their types.
        If an invalid type is recognized, throws a ValueError

        Parameters
        ----------
        r : str or numpy.ndarray
            path to a numpy array or numpy array
        """

        if type(r) == str and r[-3:] == 'npy':
            r = np.load( r )
            return r
        elif type(r) == str and r[-3:] == 'pkl':
            pin = open( r, "rb" )
            r = pickle.load(pin)
            return r
        elif type(r) == np.ndarray:
            return r
        else:
            raise ValueError(f"# Provided {tp} rates have to be either numpy.ndarray, a pickle or a string.")
    #----------------------------------------------------------------------------------------------------------------

    def build_toyexp( self ):

        """
        Builds toy experimental data based on simulated ones.
        """

        tmp       = np.copy( self.rex )
        self.rex  = np.mean( tmp, axis = -1 )
        self.eex  = np.std(  tmp, axis = 2 ) / np.sqrt( tmp.shape[2] )

    #----------------------------------------------------------------------------------------------------------------

    def load_methyl_list( self, input ):

        pin      = open(input, "rb")
        self.mnl = pickle.load(pin)

        if len(self.idx) > 0:
            mnl = []

            for j in range(self.rex.shape[-1]):
                if j not in self.idx:
                    mnl.append(self.mnl[j])

            self.mnl = mnl
    # ----------------------------------------------------------------------------------------------------------------

    def ignore_methyls( self ):

        """
        Removes from the dataset a list of specified methyls that need to be ignored.

        Parameters
        ----------
        idx : list
            indices of methyls to be removed
        """

        rex = np.empty( (self.rex.shape[0], self.rex.shape[1] - len(self.idx)) )
        rmd = np.empty( (self.rmd.shape[0], self.rmd.shape[1] - len(self.idx), self.rmd.shape[2]))
        eex = np.empty( (self.eex.shape[0], self.rex.shape[1] - len(self.idx)) )

        k = 0
        for j in range(self.rex.shape[-1]):
            if j not in self.idx:
                rex[:,k]   = self.rex[:,j]
                eex[:,k]   = self.eex[:,j]
                rmd[:,k,:] = self.rmd[:,j,:]
                k += 1

        self.rex = rex
        self.eex = eex
        self.rmd = rmd
    #----------------------------------------------------------------------------------------------------------------

    def _chi2( self, r, w ):

        """
        Computes chi squared between simulated and experimental rates.
        Reweighting is automatically accounted for.

        Parameters
        ----------
        r : int
            rate index. r=-1 means that all rates are considered in the calculation.
        w : numpy.array
            array of weights.
        """

        #rrw = np.average( self.rmd, weights = w, axis = -1 )
        rrw = np.dot( self.rmd, w[:,np.newaxis] )[:,:,0] #removes last, useless axis
        if r == -1:
            er2 = self.eex**2 + self.emd**2
            return np.sum( (self.rex - rrw)**2 / er2 )
        else:
            er2 = self.eex[r]**2 + self.emd[r]**2
            return np.sum( (self.rex[r] - rrw[r])**2 / er2 )
    #----------------------------------------------------------------------------------------------------------------

    def _rmsd( self, r, w ):

        """
        Computes chi squared between simulated and experimental rates.
        Reweighting is automatically accounted for.

        Parameters
        ----------
        r : int
            rate index. r=-1 means that all rates are considered in the calculation.
        w : numpy.array
            array of weights.
        """

        #rrw = np.average( self.rmd, weights = w, axis = -1 )
        rrw = np.dot( self.rmd, w[:,np.newaxis] )[:,:,0] #removes last, useless axis
        if r == -1:
            #er2 = self.eex**2 + self.emd**2
            return np.sum( (self.rex - rrw)**2 )
        else:
            #er2 = self.eex[r]**2 + self.emd[r]**2
            return np.sum( (self.rex[r] - rrw[r])**2 )
    #----------------------------------------------------------------------------------------------------------------

    def chi2r( self, r, w ):

        """
        Computes reduced chi squared between simulated and experimental rates.
        Reweighting is automatically accounted for.

        Parameters
        ----------
        r : int
            rate index. r=-1 means that all rates are considered in the calculation.
        w : numpy.array
            array of weights.
        """

        if r == -1:
            x2 = self._chi2( r, w )
            return x2 / ( self.m * self.r )  # divide it by the number of methyls times the number of rates
        else:
            x2 = self._chi2( r, w )
            return x2 / self.m               # divide it by the number of methyls
    #----------------------------------------------------------------------------------------------------------------

    def phi_eff( self, w ):

        """
        Computes the fraction of effective parameters.

        Parameters
        ----------
        w : numpy.ndarray
            array of weights to be optimized

        Returns
        -------
        phi : numpy.ndarray
            array of effective fraction of frames
        """

        idxs = np.where( w > 1e-50 )
        srel = np.sum( w[idxs] * np.log( w[idxs] / self.w0[idxs] ) )
        phi  = np.exp( -srel )
        return phi
    #----------------------------------------------------------------------------------------------------------------

    def _penalty( self, w, r, theta ):

        """
        Computes the value of the penalty function for a single rate.

        Parameters
        ----------
        w : numpy.ndarray
            array of weights to be optimized
        r : int
            rate index.
        theta : int
            fudge parameter.

        Returns
        -------
        r : float
            value of the penalty function given the set of weights
        """

        w   /= np.sum( w )                                           # normalize weights
        idxs = np.where( w > 1e-50 )                                 # isolate non-zero weights
        chi2 = self._chi2(r, w)                                      # compute chi2
        srel = np.sum( w[idxs] * np.log( w[idxs] / self.w0[idxs] ) ) # compute relative entropy
        p    = 0.5 * chi2 + theta * srel                             # sum all the contributions

        return p
    #----------------------------------------------------------------------------------------------------------------

    def _save( self ):

        """
        Saves reweighting results as a dictionary. The dictionary has thetas as keys and weights as arguments.
        """

        results = dict.fromkeys( self.ths, [] )

        k = 0
        for r in results.keys():
            results[r] = self.res[k].x / np.sum( self.res[k].x ) # convert current results into a dictionary
            k += 1

        # save all the results in a pickle
        with open( self.out + ".pkl", "wb" ) as fp:
            pickle.dump( results, fp )

        self.res = results
    #----------------------------------------------------------------------------------------------------------------

    def reweight( self, r ):

        """
        Reweights the data with respect to a single rate r.
        It saves the results in a .pkl file.

        Parameters
        ----------
        r : int
            rate index.
        """

        self.res = []
        bounds = [(0,1)] * len(self.w0)                                                  # weights have to be bound between 0 and 1
        opt    = { 'maxiter': 10000, 'maxfun': 10000000 }
        flags  = []                                                                      # stores error messages of minimizations that didn't converge

        w0 = self.w0
        for t in self.ths:                                                               # run minimization for all values of theta

            print(f"# THETA: {t}         \r", end  = "")

            rs = optimize.minimize( self._penalty, w0, args = (r,t), bounds = bounds, jac = False, options = opt, method = 'L-BFGS-B' )
            self.res.append( rs )
            w0 = rs.x / np.sum(rs.x)

            if not rs.success:                                                           # some minimizations had problems!
                flags.append( [t, rs.message] )

        if flags == []:
            print("\n# Done! All minimizations terminated successfully")
        else:
            print("\n# Done! Some minimizations terminated unsuccessfully: ")
            print(flags)

        self._save()
        print( f"# Saved {self.out}.pkl" )
    #----------------------------------------------------------------------------------------------------------------

    def run( self ):

        """
        Runs reweighting for all the available rates.
        The name of the saved files will be the provided output name + '_R<num>'
        In the case of all the rates, <num> = a.
        """

        for r in range( -1, self.r ):

            if r == -1:
                print("# All rates")
                add = "_Ra"
            else:
                add = f"_R{r+1}"
                print(f"# Rate {r}")

            self.out += add
            self.reweight(r)
            self.out = self.out[:-3]
            print()
    #----------------------------------------------------------------------------------------------------------------

    def load_results( self, input ):

        """
        Load a pickle file with results.

        Parameters
        ----------
        in : str
            path to the results file.
        """

        pin      = open( input, "rb" )
        self.res = pickle.load( pin )
    #----------------------------------------------------------------------------------------------------------------

    def phix2r( self, r ):

        """
        Computes the phi_eff vs reduced chi squared curve for a provided rate.

        Parameters
        ----------
        r : int
            rate index.
        """

        self.chi = []
        self.phi = []

        for k in self.res.keys():
            w = self.res[k]
            self.phi.append( self.phi_eff( w ) )
            x2 = self.chi2r( r, w )
            self.chi.append( x2 )
    #----------------------------------------------------------------------------------------------------------------

    def create_masks( self, mx, intervals ):
        mask = np.empty(mx)
        mask.fill(1)

        for i in intervals:
            mask[i[0]:i[1]].fill(0)

        return mask, 1-mask
    # ----------------------------------------------------------------------------------------------------------------

    def evaluate_sampling( self, window_length, r, label = None, outfig = None ):

        rmd = np.copy( self.rmd )
        rex = np.copy( self.rmd )
        nsamples = int( self.b / window_length )

        A = []
        for n in range(nsamples + 1):
            A.append([n * window_length, (n + 1) * window_length])

        del A[-1]

        chi2 = []
        chi2_std = []
        for i in range(nsamples):
            L = list( combinations(A, i) )

            chi2_av = []
            for l in L:

                if l == ():
                    continue
                elif len(l) > nsamples/2:
                    continue

                mask1, mask2 = self.create_masks( nsamples * window_length, l )

                rmd = rmd.compress( mask1, axis = -1 ) # apply the mask
                emd = np.std( rmd, axis = -1 ) / np.sqrt( rmd.shape[-1] )
                rmd = np.average( rmd, axis = -1 )  # average over blocks

                rex = rex.compress( mask2, axis = -1 ) # apply the opposite mask
                eex = np.std( rex, axis = -1 ) / np.sqrt( rex.shape[-1] )
                rex = np.average( rex, axis = -1 ) # average over blocks

                chi =  np.sqrt( np.sum( (rex[r] - rmd[r])**2 ) / self.m )
                chi2_av.append( chi )

                rmd = np.copy( self.rmd )
                rex = np.copy( self.rmd )

            if chi2_av == []:
                continue

            chi2.append( np.average(chi2_av) )
            chi2_std.append( np.std(chi2_av) / np.sqrt( len(chi2_av) ) )

        plt.figure(figsize=(9.55, 5))
        plt.plot(np.arange(1, nsamples / 2 + 1, 1), chi2, 'o-', c='tab:blue', markersize=8,
                 linewidth=2, markeredgecolor='k')
        plt.errorbar( np.arange(1, nsamples / 2 + 1, 1), chi2, yerr = chi2_std, c = 'tab:blue' )

        if label != None:
            plt.ylabel('RMSD R' + label + r' [s$^{-1}$]')
        else:
            plt.ylabel(r'RMSD [s$^{-1}$]')
        plt.xlabel(r'$N$ left out')


        if outfig != None:
            plt.tight_layout()
            plt.savefig(outfig + '.pdf', format='pdf')
        else:
            plt.show()
    # ----------------------------------------------------------------------------------------------------------------

    def plot_phix2r( self, r, outfig = None ):

        """
        Plots the phi_eff vs reduced chi squared curve for a provided rate.

        Parameters
        ----------
        r : int
            rate index.
        outfig : str
            path to the figure to save. If not provided, the figure will be prompted on screen and not saved.
        """

        self.phix2r( r )

        plt.figure( figsize = (5,3) )
        plt.scatter( self.phi, self.chi, c = 'tab:red', edgecolor = 'k', zorder = 10 )
        plt.plot( self.phi, self.chi, '-', c = 'tab:red' )

        plt.ylabel( r'$\chi^2_R$' )
        plt.xlabel(r'$\phi_{eff}$')
        plt.tight_layout()

        if outfig != None:
            plt.savefig( outfig + '.pdf', format = 'pdf' )
            print( f"# Saved {outfig}.pdf" )
        else:
            plt.show()
    #----------------------------------------------------------------------------------------------------------------

    def plot_phicorr( self, r, outfig = None ):

        """
        Plots the phi_eff vs Pearson correlation curve for a provided rate.

        Parameters
        ----------
        r : int
            rate index.
        outfig : str
            path to the figure to save. If not provided, the figure will be prompted on screen and not saved.
        """

        self.cor = []
        self.phi = []

        for k in self.res.keys():
            w = self.res[k]
            rrw = np.average( self.rmd, weights = w, axis = -1 )
            c   = np.corrcoef( rrw[r], self.rex[r] )[0,1]
            self.phi.append( self.phi_eff( w ) )
            self.cor.append( c )

        plt.figure( figsize = (8,5) )
        plt.scatter( self.phi, self.cor, c = 'tab:blue', edgecolor = 'k', zorder = 10 )
        plt.plot( self.phi, self.cor, '-', c = 'tab:blue' )

        plt.xlabel(r'$\phi_{eff}$')
        plt.ylabel('Pearson correlation')
        plt.tight_layout()

        if outfig != None:
            plt.savefig( outfig + '.pdf', format = 'pdf' )
            print( f"# Saved {outfig}.pdf" )
        else:
            plt.show()

        pass
    #----------------------------------------------------------------------------------------------------------------

    def plot_comparison( self, r, opt_theta = None, rate_label = None, outfig = None ):

        """
        Plots the comparison between experimental and simulated data, adding also reweighted data when provided.

        Parameters
        ----------
        r : int
            rate index.
        opt_theta : int
            theta corresponding to the optimal set of weights. If not provided, no reweighted results will be shown.
            Default: None.
        rate_label : str
            label to print on the y-axis. Is expected to be in a form similar to: r'(D$_y$)'. If not provided, a default label will be printed.
            Default: None.
        outfig : str
            path to the figure to save. If not provided, the figure will be prompted on screen and not saved.
            Default = None.
        """

        fig, ax = plt.subplots(1, 1, figsize = (7,6.5) )

        if rate_label == None:
            rate_label = r''
        chi_label   = r'$\chi^2_R=$'
        cor_label   = r'$\rho=$'
        theta_label = r'$\theta=$'
        x           = np.linspace(-1000,1000,2)
        chi0        = self.ix2[r]

        ax.scatter( self.rex[r], self.rav[r], color = 'k',       marker = 's', label = f'MD, {chi_label}{chi0:.2f}'  )
        ax.plot( x, x,'k--', zorder = -1 )

        # plot reweighting results only if optimal theta is provided
        if opt_theta != None:
            w           = self.res[opt_theta]
            opt_id      = np.where( self.ths == opt_theta )[0][0]
            rrw         = np.average( self.rmd, weights = w, axis = -1 )

            self.phix2r( r )
            chi    = self.chi2r( r, w )

            ax.scatter( self.rex[r], rrw[r],      color = 'tab:red', marker = 'o', s = 70, alpha = 0.8, \
                       label = f'AbsurdER, {theta_label}{opt_theta}, {chi_label}{chi:.2f}', edgecolor = 'k' )

            insax = ax.inset_axes([0.05,0.6,0.4,0.38])

            insax.plot( self.phi, self.chi, 'o-', c = 'tab:grey', markersize = 4, mec = 'k')
            insax.scatter( self.phi[opt_id], self.chi[opt_id], marker = 'X', c = 'tab:red', zorder = 10, s = 90, edgecolor = 'k' )
            insax.set_xlabel(r'$\phi_{eff}$', fontsize = 14)
            insax.set_ylabel(r'$\chi^2_R$', fontsize = 14)
            insax.yaxis.tick_right()
            insax.yaxis.set_label_position("right")

            insax.set_xticks([0,0.5,1])
            insax.tick_params(labelsize = 14)

        ax.set_xlabel( r'$R^{NMR}$' + rate_label + ' [s$^{-1}$]' )
        ax.set_ylabel( r'$R^{SIM}$' + rate_label + ' [s$^{-1}$]' )
        ax.set_xlim( (0, self.rex[r].max() + 5) )
        ax.set_ylim( (0, self.rex[r].max() + 5) )
        ax.legend( loc = 4 )
        plt.tight_layout()

        if outfig != None:
            plt.savefig( outfig + '.pdf', format = 'pdf' )
            print( f"# Saved {outfig}.pdf" )
        else:
            plt.show()
    #----------------------------------------------------------------------------------------------------------------

    def load_specdens( self, jex, jws, jmd ):

        """
        Loads into the class the files needed to plot the spectral densities

        Parameters
        ----------
        jex : str
            path to the pickle with the experimental spectral density functions.
        jws : int
            path to the pickle with the values of J(0), J(w) and J(2w).
        jmd : str
            path to the pickle with the simulated spectral density functions.
        """

        pin = open( jex, "rb" )
        self.jex = pickle.load( pin )
        self.jex = self.jex.T #remember to remove this with new data!!

        pin = open( jws, "rb" )
        self.jws = pickle.load( pin )
        self.jws = self.jws.T  # remember to remove this with new data!!

        pin = open( jmd, "rb" )
        self.jmd = pickle.load( pin )
        self.jmd = self.jmd.T  # remember to remove this with new data!!

        self.ignore_specdens() #ignore a set of methyl groups
        self.specdens_load = True
    #-----------------------------------------------------------------------------------------------------------------

    def ignore_specdens( self ):

        """
        Removes from the dataset of spectral densities a list of specified methyls that need to be ignored.

        Parameters
        ----------
        idx : list
            indices of methyls to be removed
        """

        jex = np.empty( (self.jex.shape[0], self.jex.shape[1] - len(self.idx), self.jex.shape[2]) )
        jmd = np.empty( (self.jmd.shape[0], self.jmd.shape[1] - len(self.idx), self.jmd.shape[2]) )
        jws = np.empty( (self.jws.shape[0], self.jws.shape[1] - len(self.idx), self.jws.shape[2]) )

        k = 0
        for j in range(self.jex.shape[1]):
            if j not in self.idx:
                jex[:,k,:] = self.jex[:,j,:]
                jmd[:,k,:] = self.jmd[:,j,:]
                jws[:,k,:] = self.jws[:,j,:]
                k += 1

        self.jex = jex
        self.jws = jws
        self.jmd = jmd
    # -----------------------------------------------------------------------------------------------------------------

    def load_rotamers( self, exrot, mdrot, ami ):

        """
        Loads into the class the files needed to plot the rotamer distributions

        Parameters
        ----------
        exrot : str
            path to the pickle with the experimental rotamer distributions.
        mdrot : int
            path to the pickle with the simulated rotamer distributions.
        ami : str
            path to the pickle with the amino acid employed in the calculation of the different rotamers.
        """

        pin = open( exrot, "rb" )
        self.exrot = pickle.load( pin )

        pin = open( mdrot, "rb" )
        self.mdrot = pickle.load( pin )

        pin = open( ami, "rb" )
        self.ami = pickle.load( pin )

        self.rot_load = True
    #-----------------------------------------------------------------------------------------------------------------

    def plot_specdens( self, idx, wd, opt_theta = None, rate_labels = [], outfig = None ):

        """
        Plots the spectral density corresponding to a specific methyl

        Parameters
        ----------
        idx : int
            methyl group index.
        wd: float
            Larmor frequency of 2H at the used magnetic field strength in MHz (Ex. 145.858415 for 2H at 950 MHz magnetic field strength)
        opt_theta : int
            theta corresponding to the optimal set of weights. If not provided, no reweighted results will be shown.
            Default: None.
        methyl_name: str
            name of the methyl group, used as a figure title
        rate_labels : list
            list of labels to print on the x-axis. Is expected to be in a form similar to: [r'(D$_y$)]'. If not provided, a default label will be printed.
            Default: [].
        outfig : str
            path to the figure to save. If not provided, the figure will be prompted on screen and not saved.
            Default = None.
        """

        if not self.specdens_load:
            raise ValueError("Spectral densities have not been loaded. Use load_specdens() for that.")

        pico    = 1. * 10 ** (-12)  # picoseconds
        omega_D = 2. * np.pi * wd * 1000000
        b       = 14  # gridpsec dimension
        wd      = np.array([0, omega_D, 2 * omega_D])
        freq    = np.linspace(0, 2 * 10 ** 9, 100)

        jws = np.average( self.jws, axis = -1 ) / pico
        jex = np.average( self.jex, axis = -1 ) / pico
        jmd = np.average( self.jmd, axis = -1 ) / pico

        fig = plt.figure( figsize = [12, 5], constrained_layout = True)
        gs  = fig.add_gridspec(1, b + self.r)

        ax1 = fig.add_subplot( gs[0, 0:b] )
        ax1.plot( wd,   jws[:, idx],       c = 'k',        ls ='',              zorder=11, marker='v', markersize=10 )
        ax1.plot( freq, jmd[:, idx], lw=4, c = 'tab:grey', ls =':', label='MD', zorder=10)
        ax1.plot( freq, jex[:, idx], lw=3, c = 'k',                 label='NMR')

        if opt_theta != None:
            jrw = np.average( self.jmd, weights = self.res[opt_theta], axis = -1 ) / pico
            rrw = np.average( self.rmd, weights = self.res[opt_theta], axis = -1 )
            ax1.plot( freq, jrw[:, idx], lw = 4, c = 'tab:red', ls = '--', label = 'ABSURDer' )

        ax1.set_ylabel(r'$J$ [ps]')
        ax1.set_xlabel(r'$\omega $ [s$^{-1}$]' + f'\n\n')
        ax1.set_yscale('log')
        ax1.xaxis.offsetText.set_fontsize(14)
        plt.legend( loc = 'upper right' )

        i = b
        j = b + 1
        for r in range(self.r):

            if rate_labels == []:
                rate = 'Rate'
            else:
                rate = 'R' + rate_labels[r]

            ax = fig.add_subplot( gs[0, i:j] )
            ax.errorbar( [rate], self.rex[r, idx], yerr = self.eex[r, idx], elinewidth = 1.2, capthick = 1.2, capsize = 3, marker = 'D',
                         markersize = 10, markeredgecolor = 'k', color = 'k')
            ax.errorbar( [rate], self.rav[r, idx], yerr = self.emd[r, idx], elinewidth = 1.2, capthick = 1.2, capsize = 3, ecolor = 'k', marker = 's',
                        markersize = 10, markeredgecolor = 'k', color = 'tab:grey')

            if opt_theta != None:
                ax.errorbar( [rate], rrw[r, idx], yerr = self.emd[r, idx], elinewidth = 1.2, capthick = 1.2, capsize = 3, ecolor = 'k', marker = 'o',
                            markersize = 11, markeredgecolor = 'k', color = 'tab:red')
            if r == 0:
                ax.set_ylabel(r'Relaxation rate [s$^{-1}$]')

            i += 1
            j += 1

        plt.suptitle(self.mnl[idx], fontsize = 20)

        if outfig != None:
            plt.savefig(outfig + '.pdf', format='pdf')
            print(f"# Saved {outfig}.pdf")
    #-----------------------------------------------------------------------------------------------------------------

    def plot_rate_distributions( self, idx, opt_theta = None, rate_labels = [], outfig = None ):

        """
        Plots the rate distributions over the blocks for a given methyl group.

        Parameters
        ----------
        idx : int
            methyl group index.
        opt_theta : int
            theta corresponding to the optimal set of weights. If not provided, no reweighted results will be shown.
            Default: None.
        methyl_name: str
            name of the methyl group, used as a figure title
        rate_labels : list
            list of labels to print on the x-axis. Is expected to be in a form similar to: [r'(D$_y$)]'. If not provided, a default label will be printed.
            Default: [].
        outfig : str
            path to the figure to save. If not provided, the figure will be prompted on screen and not saved.
            Default = None.
        """

        mpl.rcParams['xtick.labelsize'] = 18
        mpl.rcParams['ytick.labelsize'] = 18
        mpl.rcParams['axes.labelsize'] = 18
        mpl.rcParams['legend.fontsize'] = 18

        #fig = plt.figure( figsize = (20, 6) )
        fig = plt.figure(figsize=(9.55, 3))

        for r in range( self.r ):
            x      = np.linspace( self.rmd[r, idx, :].min(), self.rmd[r, idx, :].max(), num = 100 )
            kde_md = scs.gaussian_kde( self.rmd[r, idx, :], bw_method = "silverman" )
            kde_md.set_bandwidth( kde_md.scotts_factor() / 1.5 )
            kde_md = kde_md.evaluate( x )

            if opt_theta != None:
                kde_rw = scs.gaussian_kde( self.rmd[r, idx, :], bw_method = "silverman", weights = self.res[opt_theta] )
                kde_rw.set_bandwidth( kde_rw.scotts_factor() / 1.5 )
                kde_rw = kde_rw.evaluate( x )
                rrw = np.average( self.rmd, weights = self.res[opt_theta], axis = -1 )

                my = max([max(kde_md), max(kde_rw)])
            else:
                my = max(kde_md)

            my += 0.05 * my

            plt.subplot(1, self.r, r + 1)

            plt.plot( x, kde_md, lw = 3, color = 'tab:grey', zorder = -10 )
            plt.fill_between( x, kde_md, color = 'tab:grey', alpha = 0.3, label = 'MD' )
            plt.vlines( self.rav[r, idx], 0, my, color = 'tab:grey', lw = 4, linestyle = ':', zorder = 10, label = 'Average MD' )
            plt.vlines( self.rex[r, idx], 0, my, lw = 3, zorder = 1, label = 'NMR' )
            plt.axvspan( self.rex[r, idx] - self.eex[r, idx], self.rex[r, idx] + self.eex[r, idx], 0.05, 0.96, color = 'k', alpha = 0.4, zorder = 0 )

            if opt_theta != None:
                plt.plot( x, kde_rw, lw = 3, color = 'tab:red', zorder = -10)
                plt.fill_between( x, kde_rw, color = 'tab:red', alpha = 0.3, label = 'ABSURDer' )
                plt.vlines( rrw[r, idx], 0, my, color = 'tab:red', lw = 4, linestyle = '--', zorder = 5, label = 'Average ABSURDer' )

            if rate_labels == []:
                label = r'Rate [s$^{-1}$]'
            else:
                label = 'R' + rate_labels[r] + r' [s$^{-1}$]'

            plt.xlabel( label )
            if r == 0:
                plt.ylabel( 'p(R)' )
            elif r == self.r - 1:
                #plt.legend( bbox_to_anchor=(1.05, 1), loc = 'upper left', borderaxespad = 0. )
                pass

            plt.suptitle( self.mnl[idx], fontsize = 18 )
            plt.tight_layout(rect=[0, 0, 1, 0.95])

            if outfig != None:
                plt.savefig(outfig + '.pdf', format='pdf')
                print(f"# Saved {outfig}.pdf")
    #-----------------------------------------------------------------------------------------------------------------

    def plot_rotamer_distributions( self, idx, nblocks, block_size, ntrajs, opt_theta = None, outfig = None ):

        """
        Plots the rotamer distributions for a given methyl group.

        Parameters
        ----------
        idx : str
            residue name and number (ex. ILE9).
        nblocks : int
            number of blocks employed in the calculation.
        block_size : int
            size of blocks in ps.
        ntrajs : int
            number of trajectories used to compute the rotamers.
        opt_theta : int
            theta corresponding to the optimal set of weights. If not provided, no reweighted results will be shown.
            Default: None.
        outfig : str
            path to the figure to save. If not provided, the figure will be prompted on screen and not saved.
            Default = None.
        """

        if not self.rot_load:
            raise ValueError("Rotamers have not been loaded. Use load_rotamers() for that.")

        def get_hist(nblocks, blocksize, ang_methyls, mn=-180, mx=180):

            histograms = []
            for b in range(nblocks - 1):
                out = ang_methyls[b * block_size + 1:(b + 1) * block_size]
                h, _ = np.histogram(out, bins=100, range=(mn, mx))
                histograms.append(h)

            return histograms

        chi1      = ['ILE', 'LEU', 'MET', 'THR']
        chi2      = ['ILE', 'LEU', 'MET']
        ang_names = [r"$\chi_1$", r"$\chi_2$", "$\phi$", "$\psi$"]
        rng_max   = [240, 240, 180, 180]
        rng_min   = [-120, -120, -180, -180]
        shift     = [17, 17, 0, 0]
        len_traj  = int(nblocks / ntrajs)

        if idx[:3] in chi2:
            a = 2
            b = 2
            size = (15, 10)
        elif idx[:3] in chi1 and idx[:3] not in chi2:
            a = 1
            b = 3
            size = (15, 5)
        else:
            a = 1
            b = 2
            size = (12, 5)

        fig, axs = plt.subplots(a, b, figsize=size)
        angs = []
        for ax in fig.axes:
            for angg in range(4):
                if idx in self.ami[angg] and angg not in angs:
                    ang = angg
                    angs.append(ang)
                    break

            ind         = self.ami[ang].index(idx)
            tmp_exp     = self.exrot[ang][:, ind]
            hist_exp, _ = np.histogram(tmp_exp, bins=100, range=(-180, 180))
            norm        = np.sum(hist_exp)
            hist_exp    = hist_exp / norm / 3.6
            hist_exp    = np.roll(hist_exp, shift[ang])  # shifts histogram to optimal range

            tmp_md = self.mdrot[ang][:, ind]
            hist   = get_hist( nblocks, block_size, tmp_md )

            conc = ()
            for n in range(1, ntrajs + 1):
                conc = conc + (hist[(n - 1) * len_traj:n * len_traj - 1],)
            hist = np.concatenate(conc)

            hist_sum = np.average(hist, axis=0) * len(hist)
            norm     = np.sum(hist_sum)
            hist_md = hist_sum / norm / 3.6
            hist_md = np.roll(hist_md, shift[ang])

            if opt_theta != None:
                hist_sum = np.average(hist, axis=0, weights=self.res[opt_theta]) * len(hist)
                norm     = np.sum(hist_sum)
                hist_rw  = hist_sum / norm / 3.6
                hist_rw  = np.roll(hist_rw, shift[ang])

            ax.plot(np.linspace(rng_min[ang], rng_max[ang], 100), hist_exp, c='k', lw=4, label='NMR')
            ax.plot(np.linspace(rng_min[ang], rng_max[ang], 100), hist_md, c='tab:grey', lw=4, ls=':', label='MD')

            if opt_theta != None:
                ax.plot(np.linspace(rng_min[ang], rng_max[ang], 100), hist_rw, c='tab:red', lw=4, ls='--', label='ABSURDer')
            ax.set_xlabel(ang_names[ang])
            ax.set_ylabel(None)
            ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:1.2f}'))

        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        fig.suptitle(idx, fontsize=22)

        if outfig != None:
            plt.savefig(outfig + '.pdf', format='pdf')
            print(f"# Saved {outfig}.pdf")

    #------------------------------------------------------------------------------------------------------------------

    def plot_single_rotamer( self, idx, ang, nblocks, block_size, ntrajs, opt_theta = None, outfig = None ):

        """
        Plots the rotamer distributions for a given methyl group.

        Parameters
        ----------
        idx : str
            residue name and number (ex. ILE9).
        nblocks : int
            number of blocks employed in the calculation.
        block_size : int
            size of blocks in ps.
        ntrajs : int
            number of trajectories used to compute the rotamers.
        opt_theta : int
            theta corresponding to the optimal set of weights. If not provided, no reweighted results will be shown.
            Default: None.
        outfig : str
            path to the figure to save. If not provided, the figure will be prompted on screen and not saved.
            Default = None.
        """

        def get_hist(nblocks, blocksize, ang_methyls, mn=-180, mx=180):

            histograms = []
            for b in range(nblocks - 1):
                out = ang_methyls[b * block_size + 1:(b + 1) * block_size]
                h, _ = np.histogram(out, bins=100, range=(mn, mx))
                histograms.append(h)

            return histograms

        chi1      = ['ILE', 'LEU', 'MET', 'THR']
        chi2      = ['ILE', 'LEU', 'MET']
        ang_names = [r"$\chi_1$", r"$\chi_2$", "$\phi$", "$\psi$"]
        rng_max   = [240, 240, 180, 180]
        rng_min   = [-120, -120, -180, -180]
        shift     = [17, 17, 0, 0]
        len_traj  = int(nblocks / ntrajs)

        plt.figure( figsize = ( 9.55, 5 ) )
        plt.title( idx, fontsize=18, weight = 'bold')
        angs = []

        ind         = self.ami[ang].index(idx)
        tmp_exp     = self.exrot[ang][:, ind]
        hist_exp, _ = np.histogram(tmp_exp, bins=100, range=(-180, 180))
        norm        = np.sum(hist_exp)
        hist_exp    = hist_exp / norm / 3.6
        hist_exp    = np.roll(hist_exp, shift[ang])  # shifts histogram to optimal range

        tmp_md = self.mdrot[ang][:, ind]
        hist   = get_hist( nblocks, block_size, tmp_md )

        conc = ()
        for n in range(1, ntrajs + 1):
            conc = conc + (hist[(n - 1) * len_traj:n * len_traj - 1],)
        hist = np.concatenate(conc)

        hist_sum = np.average(hist, axis=0) * len(hist)
        norm     = np.sum(hist_sum)
        hist_md = hist_sum / norm / 3.6
        hist_md = np.roll(hist_md, shift[ang])

        if opt_theta != None:
            hist_sum = np.average(hist, axis=0, weights=self.res[opt_theta]) * len(hist)
            norm     = np.sum(hist_sum)
            hist_rw  = hist_sum / norm / 3.6
            hist_rw  = np.roll(hist_rw, shift[ang])

        plt.plot(np.linspace(rng_min[ang], rng_max[ang], 100), hist_exp, c='k', lw=4, label='NMR')
        plt.plot(np.linspace(rng_min[ang], rng_max[ang], 100), hist_md, c='tab:grey', lw=4, ls=':', label='MD')

        if opt_theta != None:
            plt.plot(np.linspace(rng_min[ang], rng_max[ang], 100), hist_rw, c='tab:red', lw=4, ls='--', label='ABSURDer')
        plt.xlabel(ang_names[ang])
        plt.ylabel(r'$p($' + ang_names[ang] + r'$)$')
        #ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:1.2f}'))

        if outfig != None:
            plt.tight_layout()
            plt.savefig(outfig + '.pdf', format='pdf')
            print(f"# Saved {outfig}.pdf")

    #------------------------------------------------------------------------------------------------------------------

    def plot_2d_rotamers( self, idx, nblocks, block_size, ntrajs, opt_theta, outfig = None ):

        """
        Plots the chi1-chi2 rotamer distribution for a given methyl group.

        Parameters
        ----------
        idx : str
            residue name and number (ex. ILE9).
        nblocks: int
            number of blocks employed in the calculation.
        block_size : int
            size of blocks in ps.
        opt_theta : int
            theta corresponding to the optimal set of weights. If not provided, no reweighted results will be shown.
        outfig : str
            path to the figure to save. If not provided, the figure will be prompted on screen and not saved.
            Default = None.
        """

        if not self.rot_load:
            raise ValueError("Rotamers have not been loaded. Use load_rotamers() for that.")

        ind         = self.ami[0].index(idx)
        ind2        = self.ami[1].index(idx)
        tmp_md_chi1 = self.mdrot[0][:, ind]
        tmp_md_chi2 = self.mdrot[1][:, ind2]
        len_traj    = int(nblocks / ntrajs)

        hist = []
        for b in range(nblocks - 1):
            out_chi1 = tmp_md_chi1[b * block_size + 1:(b + 1) * block_size]
            out_chi2 = tmp_md_chi2[b * block_size + 1:(b + 1) * block_size]
            h, _, _ = np.histogram2d(out_chi1, out_chi2, bins=100, range=[[-180, 180], [-180, 180]])
            hist.append(h)

        conc = ()
        for n in range(1, ntrajs + 1):
            conc = conc + (hist[(n - 1) * len_traj:n * len_traj - 1],)
        hist = np.concatenate(conc)

        hist_md = np.average(hist, axis=0) * len(hist)
        hist_md = np.roll(hist_md, 17, axis=0)
        hist_md = np.roll(hist_md, 17, axis=1)

        norm = np.sum(hist_md)
        hist_md = hist_md / norm / 3.6 / 3.6

        chi1_exp = self.exrot[0][:, ind]
        chi2_exp = self.exrot[1][:, ind2]

        h, _, _ = np.histogram2d(chi1_exp, chi2_exp, bins=100, range=[[-180, 180], [-180, 180]])
        norm = np.sum(h)
        h = h / norm / 3.6 / 3.6
        h = np.roll(h, 17, axis=0)
        h = np.roll(h, 17, axis=1)

        hist_rw = np.average(hist, axis=0, weights=self.res[opt_theta]) * len(hist)
        hist_rw = np.roll(hist_rw, 17, axis=0)
        hist_rw = np.roll(hist_rw, 17, axis=1)

        norm = np.sum(hist_rw)
        hist_rw = hist_rw / norm / 3.6 / 3.6

        plt.figure(figsize=(19, 6))

        oldcmp = mpl.cm.get_cmap('Reds', 512)
        newcmp = ListedColormap(oldcmp(np.linspace(0, 0.75, 384)))

        plt.subplot(1, 3, 1)
        plt.contourf(hist_md.T, 50, cmap=newcmp, zorder=1, origin='lower', extent=(-120, 240, -120, 240),
                     vmax=8e-4)
        plt.contour(hist_md.T, levels=np.arange(0, 8e-4, 1e-4), colors='k', linewidths=0.6, zorder=10, origin='lower',
                    extent=(-120, 240, -120, 240))
        plt.title('MD', fontsize=18)
        plt.xlabel(r'$\chi_1$ [deg]')
        plt.ylabel(r'$\chi_2$ [deg]')

        plt.subplot(1, 3, 2)
        plt.contourf(hist_rw.T, 50, cmap=newcmp, zorder=1, origin='lower',
                     extent=(-120, 240, -120, 240), vmax=8e-4)
        plt.contour(hist_rw.T, levels=np.arange(0, 8e-4, 1e-4), colors='k', linewidths=0.6, zorder=10, origin='lower',
                    extent=(-120, 240, -120, 240))
        plt.title('ABSURDer', fontsize=18)
        plt.xlabel(r'$\chi_1$ [deg]')

        plt.subplot(1, 3, 3)
        plt.contourf(h.T, 50, cmap=newcmp, zorder=1, origin='lower', extent=(-120, 240, -120, 240), vmax=8e-4)

        m = mpl.cm.ScalarMappable(cmap=newcmp)
        m.set_array(h)
        m.set_clim(0., 8e-4)
        cbar = plt.colorbar(m, boundaries=np.linspace(0, 8e-4, 100), ticks=[0, 2e-4, 4e-4, 6e-4, 8e-4])
        cbar.ax.set_xticklabels(
            [0, r'$2\times 10^{-4}$', r'$4\times 10^{-4}$', r'$6\times 10^{-4}$', r'$8 \times 10^{-4}$'])
        cbar.set_label('Probability Density')

        plt.contour(h.T, levels=np.arange(0, 8e-4, 1e-4), colors='k', linewidths=0.6, zorder=10, origin='lower',
                    extent=(-120, 240, -120, 240))
        plt.title('NMR', fontsize=18)
        plt.xlabel(r'$\chi_1$ [deg]')

        plt.suptitle(idx, fontsize=18, weight='bold')
        plt.tight_layout(rect=[0, 0, 1, 0.95])
    #------------------------------------------------------------------------------------------------------------------

    def phi_psi_rmsd( self, ang, nblocks, block_size, ntrajs, opt_theta ):

        def get_hist(nblocks, blocksize, ang_methyls, mn=-180, mx=180):

            histograms = []
            for b in range(nblocks - 1):
                out = ang_methyls[b * block_size + 1:(b + 1) * block_size]
                h, _ = np.histogram(out, bins = 100, range=(mn, mx))
                histograms.append(h)

            return histograms

        len_traj = int(nblocks / ntrajs)
        rmsd = []

        for res in self.ami[ang]:
            ind = self.ami[ang].index(res)
            tmp_exp = self.exrot[ang][:, ind]
            hist_exp, _ = np.histogram(tmp_exp, bins=100, range=(-180, 180))
            norm = np.sum(hist_exp)
            hist_exp = hist_exp / norm / 3.6

            tmp_md = self.mdrot[ang][:, ind]
            hist = get_hist(nblocks, block_size, tmp_md)

            conc = ()
            for n in range(1, ntrajs + 1):
                conc = conc + (hist[(n - 1) * len_traj:n * len_traj - 1],)
            hist = np.concatenate(conc)

            hist_sum = np.average(hist, axis=0) * len(hist)
            norm = np.sum(hist_sum)
            hist_md = hist_sum / norm / 3.6

            hist_sum = np.average(hist, axis=0, weights = self.res[opt_theta]) * len(hist)
            norm = np.sum(hist_sum)
            hist_rw = hist_sum / norm / 3.6

            rmsd_md = self.rmsd( hist_exp, hist_md )
            rmsd_rw = self.rmsd( hist_exp, hist_rw )
            rmsd.append( (rmsd_md - rmsd_rw) )

        return rmsd
    #------------------------------------------------------------------------------------------------------------------

    def rmsd( self, exp, md ):
        rmsd = np.sqrt( 1 / len(exp) * np.sum( ( exp - md )**2 ) )
        return rmsd
    #------------------------------------------------------------------------------------------------------------------

    def plot_delta_rmsds( self, ang, delta, label, outfig = None ):

        mpl.rcParams['xtick.labelsize'] = 18
        mpl.rcParams['ytick.labelsize'] = 18

        palette = []
        for r in self.ami[ang]:
            if 'ALA' in r:
                palette.append('tab:red')
            elif 'ILE' in r:
                palette.append('tab:brown')
            elif 'LEU' in r:
                palette.append('tab:green')
            elif 'THR' in r:
                palette.append('tab:orange')
            elif 'VAL' in r:
                palette.append('tab:blue')
            elif 'MET' in r:
                palette.append('tab:purple')

        custom_lines = [Patch(edgecolor='k', facecolor='tab:red'),
                        Patch(edgecolor='k', facecolor='tab:brown'),
                        Patch(edgecolor='k', facecolor='tab:green'),
                        Patch(edgecolor='k', facecolor='tab:orange'),
                        Patch(edgecolor='k', facecolor='tab:blue'),
                        Patch(edgecolor='k', facecolor='tab:purple')]

        labels = ['ALA', 'ILE', 'LEU', 'THR', 'VAL', 'MET']

        fig = plt.figure( figsize=(9.55, 6) )

        plt.bar( np.arange(0,len(delta),1), delta, edgecolor='k', color = palette, zorder=10 )

        plt.xlabel('Residues' )
        plt.ylabel(r'$\Delta $RMSD(' + label + ')' )
        plt.tight_layout()

        if outfig != None:
            plt.savefig( outfig + '.pdf', format = 'pdf')
        else:
            plt.show()
    #------------------------------------------------------------------------------------------------------------------
