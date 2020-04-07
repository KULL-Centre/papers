# Copyright (C) 2019 Kristoffer Enoe Johansson <kristoffer.johansson@bio.ku.dk>
"""Module to fit 2D protein unfolding"""

import numpy as np
import xlrd, lmfit, copy


class Unfold2DData(list):
    """Container class for 2D protein unfolding data

    A list of dict-elements that each describe a measured sample by the following items
    T : numpy.array
        Temperature series
    DENAT : float
        Denaturant concentration of that sample
    <curve_label> : numpy.array
        Measured curves with curve_label in NanoDsfParser.__supported_curves_sheetname.keys()
    FILE : string
        The name of the file from which the sample i loaded
    SAMPLE : int
        A sample identifier used for files with more samples
    ... : ...
        Additional items allowed
    """

    def get_curve_labels(self):
        """Get a list of curve_labels present in the data"""
        # This is bad if other parsers are implemented but there is currently no way of knowing potential curve labels. A parser base class could solve this
        ndp = NanoDsfParser() 
        return(ndp.get_curve_labels(self))
        # label_set = set()
        # for sample in self:
        #     for curve_label in ndp.__supported_curves_sheetname.keys():
        #         if curve_label in sample:
        #             label_set.add(curve_label)
        # return([*label_set])
    
    def get_temp_range(self):
        """Get the maximum temperature range of the samples in the list"""
        temp_min = np.min([sample['T'][0] for sample in self])
        temp_max = np.max([sample['T'][-1] for sample in self])
        return((temp_min,temp_max))
    
    def get_denat_range(self):
        """Get denaturant range of samples in list"""
        denat = [sample['DENAT'] for sample in self]
        return((np.min(denat),np.max(denat)))
    
    def merge_samples(self, every=1):
        """Merge all samples to 1D arrays, usefull for fitting

        Returns 4 1D-arrays of identical length containing temperature, denaturant, curve label, and observed intensity. 
        All curves are concatenated. Intensity values of numpy.NaN are discarded in the returned arrays.
        
        Parameters
        ----------
        every : integer
            Only return every <every> curve point in order to reduce the amount of data

        Returns
        ----------
        (t,d,c,i) : tuple of 4 numpy.array's
            Identically shaped arrays of temperature, denaturant, curve label, and observed intensity
        """
        curve_labels = self.get_curve_labels()
        # 1D arrays to return
        temp = np.array([])
        denat = np.array([])
        label = np.array([])
        observed = np.array([])
        for curve_label in curve_labels:
            for sample in self:
                if curve_label in sample:
                    temp = np.concatenate([temp, sample['T']])
                    denat = np.concatenate([denat, np.full(len(sample['T']), sample['DENAT'])])
                    label = np.concatenate([label, np.full(len(sample['T']), curve_label)])
                    observed = np.concatenate([observed, sample[curve_label]])
        # If the observed curves contains nan's these should be removed
        i = np.where(np.isnan(observed))[0]
        if len(i) > 0:
            temp     = np.delete(temp, i)
            denat    = np.delete(denat, i)
            label    = np.delete(label, i)
            observed = np.delete(observed, i)
        if every > 1:
            i = np.arange(0,len(temp),every)
            temp = temp[i]
            denat = denat[i]
            label = label[i]
            observed = observed[i]
        return(temp, denat, label, observed)

    def norm2temp(self, temperature, curve_labels=None, norm2one=False):
        """Scale curves to match at a given temperature
        
        Parameters
        ----------
        temperature : float
            The temperature to normalize at
        curve_labels : list of str
            The curve to normalize. Default is all curves
        norm2one : bool
            Normalize to 1 or mean signal (default) at given temperature

        Returns
        ----------
        sample_list : Unfold2DData
            A normalized copy of the list

        """
        ret = copy.deepcopy(self)
        
        if curve_labels is None:
            curve_labels = self.get_curve_labels()
        index_list = []
        value_dict = {}
        # Find the index and intensity value of the given temperature for each curve
        for sample in ret:
            i = np.argmin(np.abs(sample['T']-temperature))
            index_list.append(i)
            if np.abs(sample['T'][i]-temperature) > 1.0:
                print("WARNING: Temperature closest to normalization temperature %.1f is %.1f in %s sample %d" %
                      (temperature, sample['T'][i], sample['FILE'], sample['SAMPLE']))
            for curve_label in curve_labels:
                if curve_label in sample:
                    if not curve_label in value_dict.keys():
                        value_dict[curve_label] = []
                    value_dict[curve_label].append(sample[curve_label][i])
        # Average signal per curve
        curve_mean = {}
        for curve_label in curve_labels:
            curve_mean[curve_label] = np.nanmean(value_dict[curve_label])
        # Normalization factors per sample and curve
        # norm_fact = []
        for si in range(len(ret)):
            dnf = {}
            for curve_label in curve_labels:
                if not curve_label in ret[si].keys():
                    continue
                f = 1.0 / ret[si][curve_label][index_list[si]]
                nf = f if norm2one else f*curve_mean[curve_label]
                ret[si][curve_label] *= nf
                dnf[curve_label] = nf
            # norm_fact.append(dnf)
        # return(norm_fact)
        return(ret)

    
class NanoDsfParser:
    """Read and slice NanoDSF data files containing a number of sheets with measured curves for a number of
    samples (capilaries).

    Class variables
    ----------
    __supported_curves_sheetname : dict
        These are the supported curves and associated sheet names. The keys are the curve labels used to retrive 
        data from this class.
    """

    # Use lower case only
    __supported_curves_sheetname = {'I330':['330nm','330nm (unfolding)'], 'I350':['350nm','350nm (unfolding)'],
                                    'I330_REFOLD':['330nm (refolding)'], 'I350_REFOLD':['350nm (refolding)']}

    def __init__(self, *datafiles):
        """Constructor

        Parameters
        ----------
        **datafiles : str, optional
             Any number of data file names to load upon construction
        """
        self.__filenames = []  # List of names per loaded file
        self.__raw_temp = []   # List of temperature nd.array per file
        self.__raw_curves = [] # List of curve dicts per file. A curve dict has a element per supported curve each
                               #   containing an nd.array of the same length as the temperature array for that file
        self.__sample_map = [] # List of sample dicts per loaded sample. A sample dict contains a pointer to the
                               #   relevant temperature and curve arrays in the raw lists.

        for datafile in datafiles:
            self.load(datafile)
        
    def __get_temp_serie(self, book, *sheetnames, verbose=1):
        """Get the temperature serie from a workbook

        The temperature serie is given in all data sheet and this function reads all, check if they 
        are identical and returns the first. 

        In the given workbook sheets, this will look for (case-insensitive substring) \'temperature\' in cell B3

        Parameters
        ----------
        book : xlrd.Book object
             Book to read
        sheetnames : str
             Any number of sheet names to read temperature serie from
        verbose : int, optional
             Level of information dumped. Zero is warnings and errors only
        """
        temp_serie = None
        for sheetname in sheetnames:
            # The temperature label sould be in column number 2 row number 3
            col = book.sheet_by_name(sheetname).col_slice(1,2)
            if not "temperature" in col[0].value.lower():
                print("ERROR: Temperature is not in column 2 row 3 (cell B3) of sheet %s" % (sheetname))
                return(None)
            else:
                del col[0]
                while col[-1].value == '':
                    del col[-1]
            if temp_serie is None:
                temp_serie = np.array([cell.value for cell in col])
            else:
                if not len(col) == len(temp_serie):
                    print("ERROR: Temperature series in sheet %s is different length from that in previous sheet(s)" % (sheetname))
                    return(None)
                if not all(np.abs([cell.value for cell in col] - temp_serie) < 1e-4):
                    print("ERROR: Temperature series in sheet %s is different from that in previous sheet(s)" % (sheetname))
                    return(None)
        return(temp_serie)
    
    def __get_intensity_serie(self, book, sheetname, capilary_ids, verbose=1):
        """Get the fluorescence intensity series from a sheet

        In the given workbook sheet, this will look for capilary indices in row 1 starting column C, and the 
        (case-insensitive substring) \'fluorescence\' in row 3 starting column C.

        Parameters
        ----------
        book : xlrd.Book object
            Book to read
        sheetname : str
            The name of the sheet from which to read a fluorescence serie per capilary
        capilary_ids : list of int
            List of capilary identifier integers to look for in the sheet top row
        verbose : int, optional
            Level of information dumped. Zero is warnings and errors only
        """
            
        curve_series = {}
        ndata = None
        sheet_capilary_ids = []
        # The capilary numbers should be in first row and start from column number 3
        for cell in book.sheet_by_name(sheetname).row_slice(0,2):
            if cell.ctype == 2: # must be XL_CELL_NUMBER (value is float)
                sheet_capilary_ids.append(int(cell.value))
            else:
                sheet_capilary_ids.append(np.nan)
        
        # Note that sample_index may not match the capilary id given in the sheet
        for capilary_id in capilary_ids:
            sample_index = np.where(capilary_id == np.array(sheet_capilary_ids))[0]
            if len(sample_index) == 0:
                print("ERROR: Cannot find capilary %d in sheet %s" % (capilary_id,sheetname))
                return(None)
            elif len(sample_index) > 1:
                print("ERROR: Multiple occurences of capilary %d in sheet %s" % (capilary_id,sheetname))
                return(None)
            else:
                sample_index = sample_index[0]
            # The fluorescence label sould be from column number 3 and in row number 3
            col = book.sheet_by_name(sheetname).col_slice(2+sample_index,2)
            if "fluorescence" in col[0].value.lower():
                del col[0]
                while col[-1].value == '':
                    del col[-1]
                curve_series[capilary_id] = np.array([cell.value for cell in col])
            else:
                print("ERROR: Capilary %d in column %d sheet %s has no fluorescence keyword in row 2" %
                          (capilary_id,2+sample_index,sheetname))
                return(None)
        return(curve_series)

    def load(self, filename, verbose=1):
        """Load a NanoDSF excel file

        Parameters
        ----------
        filename : str
             File to read
        verbose : int, optional
             Level of information dumped. Zero is warnings and errors only
        """
        this_filename_index = len(self.__filenames)
        try:
            # This is slow and the on_demans argument has, unfortunately, no effect on xlsx files
            book = xlrd.open_workbook(filename, on_demand=True)
        except:
            print("ERROR: Could not open file \'%s\'" % (filename))
            return None

        # Find sheet names depending on experiment settings
        overview_sheetname = None
        curve_sheetnames = {}
        for sheetname in book.sheet_names():
            # Sheet containing overview and denaturant concentrations
            if sheetname.strip().lower() in ["overview"]:
                if not overview_sheetname is None:
                    print("WARNING: Ignoring multiple occurences of sheet name %s in %s, using first match" % (sheetname,filename))
                else:
                    overview_sheetname = sheetname
                continue
            
            # Sheet containing measured fluorescence of unfolding
            for curve_label in self.__supported_curves_sheetname.keys():
                if sheetname.strip().lower() in self.__supported_curves_sheetname[curve_label]:
                    if curve_label in curve_sheetnames:
                        print("WARNING: Ignoring multiple occurences of sheet name %s in %s, using first match" % (sheetname,filename))
                    else:
                        curve_sheetnames[curve_label] = sheetname
                    continue

        # Check if all sheets were found
        if overview_sheetname is None:
            print("ERROR: Could not find overview sheets in %s" % (filename))
            return(None)

        # Check if any curves were found
        if len(curve_sheetnames) < 1:
            print("ERROR: Could not find all any supported curves in %s" % (filename))
            return(None)

        # Read temperature series
        temp_series = self.__get_temp_serie(book, *curve_sheetnames.values())
        if temp_series is None:
            print("ERROR: There seems to be a problem with the temperature columns in %s" % (filename))
            return(None)

        # Read capilary id's from overview sheet
        capilary_ids = []
        for cell in book.sheet_by_name(overview_sheetname).col_slice(0,1):
            if cell.ctype == 2: # must be XL_CELL_NUMBER (value is float)
                capilary_ids.append(int(cell.value))
            else:
                break
        if len(capilary_ids) < 1:
            print("ERROR: Could not find any measured capilaries in sheet %s of file %s" % (overview_sheetname,filename))
            return(None)
        elif len(set(capilary_ids)) < len(capilary_ids):
            print("ERROR: Capilary list in sheet \'%s\' of file %s is redundant" % (overview_sheetname,filename))
            return(None)
        
        # Read denaturent concentrations from overview sheet
        col_name_cells = book.sheet_by_name(overview_sheetname).row(0)
        denat_col_index = None
        for ci in range(len(col_name_cells)):
            if any( [denat_string in (col_name_cells[ci]).value.lower() for denat_string in ["guhcl","denat","denaturent"]] ):
                if not denat_col_index is None:
                    print("ERROR: Multiple columns in the overview sheet of %s match a denaturent column, at least col %d and %d" % (filename,denat_col_index,ci))
                    return(None)
                else:
                    denat_col_index = ci
        if denat_col_index is None:
            print("ERROR: Could not find a denaturant column in the overview sheet of %s" % (filename))
            return(None)
        denat_cells = book.sheet_by_name(overview_sheetname).col_slice(denat_col_index,1,1+len(capilary_ids))
        if not np.all([cell.ctype==2 for cell in denat_cells]):
            print("ERROR: Not all denaturant values seems to numbers in column %d of sheet %s file %s" % (denat_col_index,overview_sheetname,filename))
            return(None)
        denat_val = [cell.value for cell in denat_cells]
        denat = dict(zip(capilary_ids,denat_val))

        # Read intensity curves
        curves = {}
        for curve_label in curve_sheetnames.keys():
            curves[curve_label] = self.__get_intensity_serie(book, curve_sheetnames[curve_label], capilary_ids, verbose)
            if curves[curve_label] is None:
                print("ERROR: There seems to be a problem with the %s sheet in %s" % (curve_sheetnames[curve_label],filename))
                return(None)

        # Read T-slope from overview sheet col 5 row 1 (F2)
        T_slope_val = [cell.value for cell in book.sheet_by_name(overview_sheetname).col_slice(5,1,1+len(capilary_ids))]
        T_slope = dict(zip(capilary_ids,T_slope_val))
            
        # At this point everything checks so update object variables
        self.__filenames.append(filename)
        self.__raw_temp.append(temp_series)
        self.__raw_curves.append(curves)

        # update __sample_map
        for capilary_id in capilary_ids:
            sample = {'T':self.__raw_temp[this_filename_index], 'DENAT':denat[capilary_id], 'TSLOPE':T_slope[capilary_id], 'FILE':self.__filenames[this_filename_index], 'SAMPLE':capilary_id}
            for curve_label in curve_sheetnames.keys():
                sample[curve_label] = self.__raw_curves[this_filename_index][curve_label][capilary_id]
            # If I get python right, the arrays will not be copied but only pointed to as long as they stay unchanged
            self.__sample_map.append(sample)
        book.release_resources()
        del book
        return(capilary_ids)

    def get(self, indices=None, temp_min=None, temp_max=None, curve_labels=None, verbose=0):
        """ Return a sliced copy of selected data 

        Parameters
        ----------
        indices : iterator returning integers, optional
             The samples to return. Default is all samples
        temp_min, temp_max : float, optional
             Temperature range to return. Default is all range
        curve_labels : iterable of curve labels
             Curves to return. Default is all curves
        verbose : int
             Level of information dumped. Zero is warnings and errors only
        """
        if curve_labels is None:
            curve_labels = self.__supported_curves_sheetname.keys()
        else:
            for label in curve_labels:
                if not label in self.__supported_curves_sheetname.keys():
                    print("ERROR: Curve label \'%s\' not supported" % (label))
                    return(None)

        if indices is None:
            indices = range(len(self.__sample_map))
        else:
            if np.min(indices) < 0 or np.max(indices) >= len(self.__sample_map):
                print("ERROR: Data range \'%s\' is outside range 0-%d" % (str(indices),len(self.__sample_map)-1))
                return(None)
        
        ret = Unfold2DData()
        for si in indices:
            # Number of data points for this sample
            ndata = len(self.__sample_map[si]['T'])
            # Determine temp_min for this sample
            if temp_min is None:
                Ti_begin = 0
            else:
                if temp_min < self.__sample_map[si]['T'][0]:
                    print("WARNING: Tmin %.1f less than Tmin of %s capilary %d - using Tmin %.1f" %
                          (temp_min,self.__sample_map[si]['FILE'], self.__sample_map[si]['SAMPLE'],self.__sample_map[si]['T'][0]))
                    Ti_begin = 0
                elif temp_min > self.__sample_map[si]['T'][ndata-1]:
                    print("WARNING: Tmin %.1f greater than Tmax %.1f of %s capilary %d - discarding sample" %
                          (temp_min, self.__sample_map[si]['T'][ndata-1], self.__sample_map[si]['FILE'], self.__sample_map[si]['SAMPLE']))
                    continue
                else:
                    Ti_begin = np.argmin(np.abs(self.__sample_map[si]['T'] - temp_min))
            this_temp_min = self.__sample_map[si]['T'][Ti_begin]
            # Determine temp_max for this sample
            if temp_max is None:
                Ti_end = ndata
            else:
                if temp_max < self.__sample_map[si]['T'][0]:
                    print("WARNING: Tmax %.1f is less than Tmin %.1f of %s capilary %d, using %.1f - discarding sample" %
                          (temp_max, self.__sample_map[si]['T'][0], self.__sample_map[si]['FILE'], self.__sample_map[si]['SAMPLE']))
                    continue
                elif temp_max > self.__sample_map[si]['T'][ndata-1]:
                    print("WARNING: Tmax %.1f is greater than Tmax of %s capilary %d - using Tmax %.1f" %
                          (temp_max, self.__sample_map[si]['FILE'], self.__sample_map[si]['SAMPLE'], self.__sample_map[si]['T'][ndata-1]))
                    Ti_end = ndata
                else:
                    Ti_end = np.argmin(np.abs(self.__sample_map[si]['T'] - temp_max))
                # if next element is equal or only slightly larger include it
                if ndata < (ndata-1):
                    if self.__sample_map[si]['T'][Ti_end]-temp_max < 1e-12:
                        Ti_end += 1
            this_temp_max = self.__sample_map[si]['T'][Ti_end-1]
            if Ti_begin >= Ti_end:
                print("ERROR: Bad temperature range [%.1f,%.1f] index [%d,%d] requested for  %s capilary %d" %
                      (this_temp_min, this_temp_max, Ti_begin, Ti_end, self.__sample_map[si]['FILE'], self.__sample_map[si]['SAMPLE']))
                return(None)
            if verbose > 0:
                print("Temperature range [%.1f,%.1f] index [%d,%d] returned for sample index %d" %
                      (this_temp_min, this_temp_max, Ti_begin, Ti_end, si))
            
            # Determine curves for this sample
            this_curve_labels = []
            for curve_label in curve_labels:
                if curve_label in self.__sample_map[si].keys():
                    this_curve_labels.append(curve_label)
                elif verbose > 0:
                    print("No curve with label \'%s\' in %s sample %d - skip curve in this sample" %
                          (curve_label, self.__sample_map[si]['FILE'], self.__sample_map[si]['SAMPLE']))
            if len(this_curve_labels) < 1:
                print("No curves left in %s sample %d - skip sample" % (self.__sample_map[si]['FILE'], self.__sample_map[si]['SAMPLE']))
                continue
            this_sample = copy.deepcopy(self.__sample_map[si])
            this_sample['T'] = this_sample['T'][Ti_begin:Ti_end]
            for curve_label in self.__supported_curves_sheetname.keys():
                if curve_label in this_curve_labels:
                    this_sample[curve_label] = this_sample[curve_label][Ti_begin:Ti_end]
                elif curve_label in this_sample:
                    del this_sample[curve_label]
            ret.append(this_sample)
        return(ret)

    def get_curve_labels(self, sample_list=None):
        """Get list of curve_labels used in a list of samples or in all loaded data

        Parameters
        ----------
        sample_list : Unfold2DData
        """
        if sample_list is None:
            sample_list = self.__sample_map
        label_set = set()
        for sample in sample_list:
            for curve_label in self.__supported_curves_sheetname.keys():
                if curve_label in sample:
                    label_set.add(curve_label)
        return([*label_set])
    
    def get_temp_range(self, sample_list=None):
        """Get maximum temperature range of a sample list or of all loaded data"""
        temp_min = np.min([sample['T'][0] for sample in self.__sample_map])
        temp_max = np.max([sample['T'][-1] for sample in self.__sample_map])
        return((temp_min,temp_max))
    
    def info(self):
        """Dump info on current content"""
        for fi in range(len(self.__filenames)):
            curve_labels = [*self.__raw_curves[fi].keys()]
            print("File %2d with %2d samples, %1d curves (%s) and %4d temperature points : \'%s\'" %
                  (fi, len(self.__raw_curves[fi][curve_labels[0]]), len(self.__raw_curves[fi]), str(curve_labels),
                   len(self.__raw_temp[fi]), self.__filenames[fi]))
        return(len(self.__sample_map))

    
class Unfold2DModel(lmfit.Model):
    """A 2-dimensional unfolding model that considers denaturant and temperature unfolding simultaneously.

    A thermodynamic function for unfolding is required upon construction and this class mostly provides
    functions for guessing parameters. All temperature are assumed to be in celsius.

    The provided function should have the form

        func(temp_target, temp, denat, label, parameter1, parameter2, ...)

    where the variables are identified in independent_vars

        Global2DUnfolding.Unfold2DModel(func, independent_vars=['variable1','variable2', ...])

    The parameter names should match the names in the guessing functions. More variables may be included 
    but then parent_fit should be used instead if fit.
    """

    def __init__(self, func, independent_vars, nan_policy='raise', **kwargs):
        """Constructor
        Parameters
        ----------
        func : fitting function
            Samples to fit. 
        params : lmfit.Parameters() instance
            Initial parameters
        temp_target : float
            Target temperature to center the parametirzed curve
        every : integer
            See Unfold2DData.merge_samples for details.
        kwargs : additional keyword arguents
            Passed to Model.fit

        See also documentation for lmfit.Model
        """
        for var in ['temp_target', 'temp', 'denat', 'label']:
            if not var in independent_vars:
                print("WARNING: Fit function does not have an independent variable called \'%s\' - use parent_fit instead of fit" % (var))
        kwargs.update({'nan_policy': nan_policy})
        super().__init__(func, independent_vars, **kwargs)        

    def parent_fit(self, *args, **kwargs):
        """New definition of the parrent Model.fit because of overwriting"""
        return(super().fit(*args, **kwargs))

    def fit(self, sample_list, params, temp_target, every=1, **kwargs):
        """Wrapper to Model.fit that flattens a Unfolding2D-style sample list into individual numpy.array's.
        
        See Unfold2DData.merge_samples for details.

        Parameters
        ----------
        sample_list : Unfold2DData
            Samples to fit. 
        params : lmfit.Parameters() instance
            Initial parameters
        temp_target : float
            Target temperature to center the parametirzed curve
        every : integer
            See Unfold2DData.merge_samples for details.
        kwargs : additional keyword arguents
            Passed to Model.fit
        
        Returns
        ----------
        result :  lmfit.ModelResult object
            A result object containing the fit
        """
        (temp,denat,label,observed) = sample_list.merge_samples(every=every)
        return(self.parent_fit(observed, params, temp_target=temp_target, temp=temp, denat=denat, label=label, **kwargs))

    def has_param(self, *param_names):
        """Test if a model has parameters of given names"""
        valid_param_names = self.make_params()
        for param_name in param_names:
            if not param_name in valid_param_names:
                return(False)
        return(True)
    
    def guess_Cp_m(self, N_res, denat="GuHCl", R=None, params=None):
        """Guess heat capacity change upon unfolding and m-value
        
        Based on the number of residues in the protein according to Geierhaas et al 2006. 
        Units are kcal/mol/K and kcal/mol/M but adapted to the provided ideal gas constant R

        Parameters
        ----------
        N_res : int
            Number of residues in the protein
        denat : string
            The denaturant in ["GuHCl","Urea"]
        R : float
            Ideal gas constant to provide energy units, Default kcal/mol
        params : lmfit.Parameters() instance, optional
            Parameter object to have dCp and m0 updated. If None, a new param object is made
        
        Returns
        ----------
        params : lmfit.Parameters() instance
            The guessed parameters
        """        
        dCp = (-130+16*N_res)/1000.0
        if denat.lower() in ["guhcl","gdmcl","gdnhcl"]:
            m = (-580+36*N_res)/1000.0
        elif denat.lower() == "urea":
            m = (70+13*N_res)/1000.0
        else:
            print("ERROR: unknown denaturent \'%s'\'" % (denat))
        # Convert units from kcal using the ideal gas constant
        R_kcal = 1.985e-3
        units = 1.0 if R is None else R/R_kcal
        dCp = dCp * units
        m = m * units
        if params is None:
            params = self.make_params()
        for (pn,pv) in zip(['dCp','m0'],[dCp,m]):
            params[pn].value = pv
        params.update_constraints()
        return params
    
    def guess_unfolded_curved(self, unfolded_sample, temp_target, curve_label, params=None):
        """Guess the second order polynomium parameters of a curved unfolded temperature-baseline

        Parameters
        ----------
        unfolded_sample : dict, Unfolding2D-style sample
            Sample that is unfolded at all temperatures from which a curved baseline is fitted
        temp_target : float
            Target temperature to center the parametirzed curve
        curve_label : string
            The label of the curve in unfolded_sample to use
        params : lmfit.Parameters() instance
            Parameter object to have a_U, b_U_T, b_U_T2 updated AND to read denaturant correction. 
            If None, a new param object is made
        
        Returns
        ----------
        params : lmfit.Parameters() instance
            The guessed parameters
        """
        # Parameter names
        denat_slope_name = 'b_U_D_'+curve_label
        param_names = [pn+'_'+curve_label for pn in ['a_U', 'b_U_T', 'b_U_T2']]
        if not self.has_param(denat_slope_name, *param_names):
            print("ERROR: Baseplane parameters for curve \'%s\' is not in model" % (curve_label))
            return(None)
        
        if params is None:
            params = self.make_params()
 
        b_U_D = params[denat_slope_name].value
        (b_U_T2,b_U_T,a_U) = np.polyfit(unfolded_sample['T']-temp_target,
                                        unfolded_sample[curve_label]-unfolded_sample['DENAT']*b_U_D, 
                                        deg=2)
        for (pn,pv) in zip(param_names,[a_U, b_U_T, b_U_T2]):
            params[pn].value = pv
        params.update_constraints()
        return params
        
    def guess_baseplane(self, sample_list, temp_target, curve_label, 
                        dT=5.0, dD=1.0, highest_denat_unfolded=True, params=None):
        """Guess folded and unfolded base-plane parameters defined as

        F(T, D; a_F, b_F_D, b_F_T)         = a_F + b_F_D D + b_F_T (T-Ttarget)
        U(T, D; a_U, b_U_D, b_U_T, b_U_T2) = a_U + b_U_D D + b_U_T (T-Ttarget) + b_U_T2 (T-Ttarget)**2
        
        Parameters
        ----------
        sample_list : list of dict, Unfolding2D style sample
            Samples from which plane parameters are guessed
        temp_target : float
            Target temperature to center the parametirzed planes
        curve_label : string
            The label of the curve in the samples to use
        dT, dD : floats
            Temperature and denaturant steps to estimate slopes from the folded and unfolded corners
        highest_denat_unfolded : bool, optional
            If true, the unfolded temperature slope is estimated from the entire unfolded curve (ignoring dT here). 
            Useful if the unfolded temperature baseline is curved
        params : lmfit.Parameters() instance, optional
            Parameter object to have dCp and m0 updated. If None, a new param object is made
        
        Returns
        ----------
        params : lmfit.Parameters() instance
            The guessed parameters
        """
        
        # Check if the requested curve_labal has baseplane parameters in the model function
        param_names = [pn+'_'+curve_label for pn in ['a_F', 'b_F_D', 'b_F_T', 'a_U', 'b_U_D', 'b_U_T']]
        if not self.has_param(*param_names):
            print("ERROR: Baseplane parameters for curve \'%s\' is not in model" % (curve_label))
            return(None)
        
        # Extract denaturant dimension and find sample with minimum denaturant concentration
        denat = np.array([d['DENAT'] for d in sample_list])
    
        ### Folded corner of base-plane
        # (D0,T0) is the corner from which the slopes are estimated to (D0,T1) and (D1,T0)
        i_D0 = np.argmin(denat)
        i_T0 = 0
        T0 = sample_list[i_D0]['T'][i_T0]
        f_F = sample_list[i_D0][curve_label][i_T0]
        if (len(denat) > 1):
            # Find D1, this sample estimation should be different from D0
            i_D1 = i_D0; e = 1e-2; D1 = denat[i_D0]+dD-e
            while i_D1 == i_D0:
                D1 += e
                i_D1 = np.argmin(np.abs(denat-D1))
            #print(i_D0,i_D1,denat)
            # Find the temperature index of T0 value in sample D1
            i_T0_D1 = np.argmin(np.abs(sample_list[i_D1]['T'] - T0))
            if np.abs(sample_list[i_D1]['T'][i_T0_D1] - T0) > 1e-2:
                print("WARNING: The folded denaturant slope is estimated at different temperature %.2f and %.2f" % 
                      (sample_list[i_D1]['T'][i_T0_D1],T0))
            b_F_D = (sample_list[i_D1][curve_label][i_T0_D1] - f_F) / (denat[i_D1] - denat[i_D0])
        else:
            b_F_D = 0.0
        # Temperature point for slope estimation, T1
        T1 = T0+dT
        i_T1 = np.argmin(np.abs(sample_list[i_D0]['T'] - T1))
        if i_T1==i_T0:
            i_T1 += 1
        # Temperature slope in folded corner
        b_F_T = (sample_list[i_D0][curve_label][i_T1] - f_F) / (sample_list[i_D0]['T'][i_T1] - T0)
        #print(i_T1, T0, sample_list[i_D0]['T'][i_T1], sample_list[i_D0]['curve_label'][curve_label][i_T1], f_F)
        # Temperature intersect at target temperature
        a_F = f_F + b_F_T*(temp_target-T0) - b_F_D*denat[i_D0]
        ### Unfolded corner of base-plane
        # (D0,T0) is the corner towards which the slopes are estimated from (D0,T1) and (D1,T0)
        i_D0 = np.argmax(denat)
        i_T0 = len(sample_list[i_D0]['T'])-1
        T0 = sample_list[i_D0]['T'][i_T0]
        f_U = sample_list[i_D0][curve_label][i_T0]
        if (len(denat) > 1):
            # Sample for denaturant slope estimation should be different
            i_D1 = i_D0; e = 1e-2; D1 = denat[i_D0]-dD+e
            while i_D1 == i_D0:
                D1 -= e
                i_D1 = np.argmin(np.abs(denat-D1))
            #print(i_D0,i_D1,D1,denat)
            # Find temperature index of the T0 value in D1
            i_T0_D1 = np.argmin(np.abs(sample_list[i_D1]['T'] - T0))
            if np.abs(sample_list[i_D1]['T'][i_T0_D1] - T0) > 1e-2:
                print("WARNING: The unfolded denaturant slope is estimated at different temperature %.2f and %.2f" % 
                      (sample_list[i_D1]['T'][i_T0_D1], T0))
            b_U_D = (f_U - sample_list[i_D1][curve_label][i_T0_D1]) / (denat[i_D0] - denat[i_D1])
        else:
            i_D1 = i_D0
            b_U_D = 0.0
        #print(i_T0,i_T0_D1,f_U,sample_list[i_D1]['curve_label'][curve_label][i_T0_D1])
        if highest_denat_unfolded:
            # Temperature slope when denatured at all tempertures
            b_U_T = (f_U - sample_list[i_D0][curve_label][0]) / (T0 - sample_list[i_D0]['T'][0])
        else:
            # Temperature slope in folded corner
            T1 = T0-dT
            i_T1 = np.argmin(np.abs(sample_list[i_D0]['T']-T1))
            if i_T1==i_T0:
                i_T1 -= 1
            b_U_T = (f_U - sample_list[i_D0][curve_label][i_T1]) / (T0 - sample_list[i_D0]['T'][i_T1])
        #print(i_T1, T0, sample_list[i_D0]['T'][i_T1], sample_list[i_D0]['curve_label'][curve_label][i_T1], f_U)
        # Temperature intersect at target temperature
        a_U = f_U - b_U_T*(T0-temp_target) - b_U_D*(denat[i_D0]-denat[i_D1])
        if params is None:
            params = self.make_params()
        for (pn,pv) in zip(param_names,[a_F, b_F_D, b_F_T, a_U, b_U_D, b_U_T]):
            params[pn].value = pv
        params.update_constraints()
        return params

    def calc_sample(self, param, temp_target, temp, denat, curve_labels):
        """Calculate a Unfolding2D-style sample given parameters and variables

        Parameters
        ----------
        params : lmfit.Parameters() instance
            Parameter of model to evaluate
        temp_target : float
            Target temperature
        temp_serie : numpy.array
            Temperature serie for model evaluation
        denat : numpy.array
            Denaturant serie for model evaluation
        curve_labels : list of strings
            Curve labels to include in the sample, must be in model
        """
        ret = {'T':temp, 'DENAT':denat, 'FILE':"calculated", 'SAMPLE':-1}
        for curve_label in curve_labels:
            curve = self.eval(param, temp_target=temp_target, temp=temp, 
                              denat=denat, label=np.full(len(temp),curve_label))
            ret[curve_label] = curve
        return(ret)

    
# Example function to fit
def unfold2d_I330_I350(temp_target, temp, denat, label,
                       a_F_I330=1.0, b_F_D_I330=0.0, b_F_T_I330=0.0, a_U_I330=1.0, b_U_D_I330=0.0, b_U_T_I330=0.0, b_U_T2_I330=0.0,
                       a_F_I350=1.0, b_F_D_I350=0.0, b_F_T_I350=0.0, a_U_I350=1.0, b_U_D_I350=0.0, b_U_T_I350=0.0, b_U_T2_I350=0.0,
                       Tm=70.0, dH=300.0, dCp=2.0, m0=1.0, m1=0.0):
    R = 8.314e-3   # Gas constant in kJ/K/mol
    K0 = -273.15   # Zero Kelvin in Celcius
    T_K = temp - K0; Tm_K = Tm - K0; ToTm = T_K/Tm_K
    dT = temp - temp_target;  dT2 = dT*dT
    # Free energy of unfolding, dH, dCp and m should be positive
    exp_dG = np.exp( -(dH*(1-ToTm) + dCp*(T_K-Tm_K-T_K*np.log(ToTm)) - denat*(m0 + m1*dT))/(R*T_K) )
    # 330nm curve_label with folded (F) and unfolded (U) base planes
    I330 = (a_F_I330 + b_F_D_I330*denat + b_F_T_I330*dT + exp_dG*(a_U_I330 + b_U_D_I330*denat + b_U_T_I330*dT + b_U_T2_I330*dT2)) / (1+exp_dG)
    # 350nm curve_label with folded (N) and denatured (D) base planes
    I350 = (a_F_I350 + b_F_D_I350*denat + b_F_T_I350*dT + exp_dG*(a_U_I350 + b_U_D_I350*denat + b_U_T_I350*dT + b_U_T2_I350*dT2)) / (1+exp_dG)
    return(np.where(label=='I330', I330, 0.0) + np.where(label=='I350', I350, 0.0))
    

# Commandline main/example
if __name__ == "__main__":
    import os,sys,readline
    filename = input("Data file: ")
    
    if not os.path.isfile(filename):
        print("ERROR: Cannot find file %s" % (filenames))
        sys.exit(-1)
    nano_dsf = NanoDsfParser()
    capilary_list = nano_dsf.load(filename)
    capilary2index = dict(zip(capilary_list,range(len(capilary_list))))
    nano_dsf.info()

    def rlinput(prompt, prefill=''):
        readline.set_startup_hook(lambda: readline.insert_text(prefill))
        try:
            return input(prompt)
        finally:
            readline.set_startup_hook()
    capilary_include = eval(rlinput('\nSelect samples to include: ', prefill=str(capilary_list)))
    index_include = [capilary2index[cap] for cap in capilary_include]
    norm_temp = float(rlinput('Set normalization temperature (negative to skip): ', prefill="80.0"))
    (Tmin,Tmax) = nano_dsf.get_temp_range()
    Tmin = float(rlinput('Set Tmin for fit: ', prefill=str(Tmin)))
    Tmax = float(rlinput('Set Tmax for fit: ', prefill=str(Tmax)))	
    T0 = float(rlinput('Set temperature origin: ', prefill="25.0"))
    Nres = float(rlinput('Number of residues in protein (negative to skip): '))
    
    data = nano_dsf.get(index_include, Tmin, Tmax, ['I330','I350'])
    if norm_temp > 0:
       data = data.norm2temp(norm_temp, norm2one=True)

    model = Unfold2DModel(unfold2d_I330_I350, independent_vars=['temp_target','temp','denat','label'])

    param_init = model.make_params()
    model.guess_baseplane(data, T0, 'I330', params=param_init)
    model.guess_baseplane(data, T0, 'I350', params=param_init)
    if Nres > 0:
        model.guess_Cp_m(Nres, params=param_init, R=8.314e-3)

    fit = model.fit(data, param_init, T0)
    print(fit.fit_report(show_correl=False))
