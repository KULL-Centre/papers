#!/storage1/shared/software/anaconda3/bin/python3

# Copyright (C) 2020 Kristoffer Enoe Johansson <kristoffer.johansson@bio.ku.dk>

"""Module for data handling in the PRISM project

This module implments classes for parsing (PrismParser) and handling 
(PrismData derived classes) of data files.

In general, calling PrismParser.read(filename) will return a data object of 
the same derived class, e.g. a VariantData object.
See documentation of derived parser and data classes for help.

See README.md for general file format definitions.
"""

__version__ = 1.001

from Bio import Seq,SeqRecord,SeqIO,pairwise2,SubsMat
from Bio.SubsMat import MatrixInfo
import numpy as np
import pandas as pd
import yaml,csv,copy,time


class PrismFormatError(Exception):
    """Exception raised for format errors in prism data files"""

class PrismValueFail(Exception):
    """Exception raised when a threshold was not met"""


class PrismParser:
    """Read and write prism data files

    Examples
    --------
    Parsing a data file
    >>> parser = PrismData.PrismParser()
    >>> data = parser.read("prism_mave/prism_mave_001_GFP.txt")

    Parsing the YAML file header only is much faster for large files
    >>> meta_data = parser.read_header("prism_mave/prism_mave_001_GFP.txt")
    """

    def read(self, filename, verbose=0):
        """Read data and meta data from a PRISM data file

        Read the meta-data section (via read_header) and data section (via Pandas.read_csv) and 
        check format adherence and the consistancy between the two.

        Parameters
        ----------
        filename : string
            Name and path of file to read
        verbose : int, optional
            Level of output

        Returns
        -------
        PrismData derived instance
        """

        if verbose > 0:
            start_time = time.process_time()
            
        # Read meta-data and check for necessary fields
        metadata = self.read_header(filename, verbose=verbose)

        # Read data section
        dataframe = pd.read_csv(filename, delim_whitespace=True, comment='#', header=0,
                             keep_default_na=True, na_values=['Na','na'])
    
        # Determine type of data
        data = None
        if dataframe.keys()[0] == "variant":
            data = VariantData(metadata, dataframe)
        # elif dataframe.keys()[0] == "resi":
        #     data = ResidueData(metadata, dataframe)
        #     checks ...
        else:
            raise PrismFormatError("Could not determine data type of file %s from first column name %s" % (filename, dataframe.columns[0]))

        data.check(verbose=verbose)
        
        if verbose > 0:
            elapsed_time = time.process_time() - start_time
            print("Total time of parsing %s %.1f sec" % (filename,elapsed_time))
            
        return(data)

    def write(self, filename, data, comment_lines=None, no_overwrite_metadata=False, verbose=0):
        """Write a PrismData object to a file

        Will check the consistency of variants and cloumns meta-data fields.

        Parameters
        ----------
        filename : string
            Name and path of file to write
        data : VariantData
            Data set to write
        comment_lines : list of strings
            Each list element is printed as a comment line between header and data sections
        no_overwrite_metadata : bool
            Do not allow the check function to overwrite meta data. Except instead
        """

        # Replace whitespace in column names if present
        has_whitespace = [cn.find(" ") > -.5 for cn in data.metadata['columns'].keys()]
        if any(has_whitespace):
            if verbose > 0:
                print("White space in columns names replaced with _ (they are white-space separated)") 
            # Replace in data frame
            data.dataframe.columns =  [cn.strip().replace(" ","_") for cn in data.dataframe.columns]
            # Replace in meta data
            data.metadata["columns"] = {k.strip().replace(" ","_"):v for (k,v) in data.metadata['columns'].items()}
        
        # Check data before printing
        data.check(no_overwrite_metadata=no_overwrite_metadata, verbose=verbose)

        # if present, remove the index columns (aa_ref, etc.)
        write_cols = list(data.dataframe.columns)
        for col in data.index_column_names:
            if col in write_cols:
                write_cols.remove(col)
                                                    
        # Write a file
        with open(filename, "w") as fh:
            fh.write(data.header_to_string(skip_fields=['filename']))
            if comment_lines:
                for line in comment_lines:
                    fh.write("# "+line+"\n")
            fh.write("#\n")
            # data.dataframe.to_csv(fh, float_format='%8.4f', sep=' ', header=True, index=False, columns=write_cols,
            data.dataframe.to_csv(fh, sep=' ', header=True, index=False, columns=write_cols,
                                  quoting=csv.QUOTE_NONE, escapechar=' ', na_rep="    NA")

    def is_aa_one_nat(self, sequence, additional=""):
        """Test if a sequence consists of single-letter natural amino acids only (case-insensitive)
        
        Parameters
        ----------
        sequence : string
            Amino acid sequence
        additional : string, optional
            Additional characters to allow, case insensitive
        """
        for aa in sequence.upper():
            if not (aa in "ACDEFGHIKLMNPQRSTVWY" or aa in additional.upper()):
                return(False)
        return(True)
    
    def read_header(self, filename, verbose=0):
        """Read the YAML header of a data file and run check_header

        Parameters
        ----------
        filename : string
            Path and name of the file to read
        verbose : int, optional
            Level of output        
        
        Returns
        -------
        dict
            The meta data as a dictonary of one or more levels of dictionaries
        """
        header_begin = False
        header_str = ""
        with open(filename, "r") as file_handle:
            for line in file_handle:
                if line.replace(' ','')[0:4] == '#---':
                    if header_begin:
                        header = yaml.safe_load(header_str)
                        header['filename'] = filename.split('/')[-1]
                        self.check_header(header)
                        return(header)
                    else:
                        header_begin = True
                        continue
                if header_begin:
                    if not line.strip()[0] == "#":
                        raise PrismFormatError("Missing comment char inside YAML header")
                    s = line.replace('#','',1)
                    header_str += s
        if header_begin:
            raise PrismFormatError("Header never ended")
        else:
            raise PrismFormatError("No header found")

    def check_header(self, header):
        """Check a header for fields required by all data files"""
        if not 'version' in header.keys():
            raise PrismFormatError("Header has no \'version\' field")
        if not 'protein' in header.keys():
            raise PrismFormatError("Header has no \'protein\' field")
        if not 'name' in header['protein'].keys():
            raise PrismFormatError("Header has no \'protein: name\' field")
        if not 'sequence' in header['protein'].keys():
            raise PrismFormatError("Header has no \'protein: sequence\' field")
        if not 'uniprot' in header['protein'].keys():
            raise PrismFormatError("Header has no \'protein: uniprot\' field")
        if 'first_residue_number' in header['protein'].keys():
            if int(header['protein']['first_residue_number']) < 0:
                raise PrismFormatError("First residue number must be non-negative")
        if not 'columns' in header.keys():
            raise PrismFormatError("Header has no \'columns\' field")
        if 'filename' in header.keys():
            data_type = header['filename'].split('_')[1]
            if not data_type.lower() in header.keys():
                raise PrismFormatError("Header has no \'%s\' field but filename indicates this data type" % (data_type))

    def write_meta_file(self, data_list, filename, output_type="csv"):
        """Write a overview file of meta-data (headers) from a list of data sets

        Parameters
        ----------
        data_list : list of PrismData-derived objects or meta data dictionaries
            The data to output
        filename : string
            Name and path of output file
        output_type : string, optional
            One of the following
            -csv: Excel-readable CSV overview of all meta-data
            -fasta: Header sequences in FASTA format
        """
        header_list = [data.metadata if isinstance(data,PrismData) else data for data in data_list]
        if output_type == "csv":
            self.__dump_header_csv(filename, header_list)
        elif output_type == "fasta":
            self.__dump_fasta(filename, header_list)
        else:
            raise ValueError("write_meta_file argument output_type must be \'csv\' or \'fasta\'")

    def __merge_header_fields(self, header_list, mode="union"):
        """Determine a common set of header fields"""
        
        # Recursive grap of keys - because we can-can
        def update_keys(key_dic, dic):
            """Update key_dic with keys from dic recursively"""
            for key in dic.keys():
                if not key in key_dic.keys():
                    key_dic[key] = {}
                if type(dic[key]) == dict:
                    update_keys(key_dic[key],dic[key])
                    
        # Read all header keys
        common_fields = {}
        for header in header_list:
            update_keys(common_fields, header)
        # Return according to mode
        if mode.lower()=="union":
            return(common_fields)
        elif mode.lower()=="intersect":
            raise NotImplementedError("__merge_header_fields mode==intersect not implemented")
        else:
            raise ValueError("Header merge mode must be \'union\' or \'intersect\'")

    def __dump_header_csv(self, filename, header_list):
        common_header_keys = self.__merge_header_fields(header_list, mode="union")
        
        # Columns fields have different names so rename to column01, column02, etc
        ncol_list = [len(header['columns'].keys()) for header in header_list]
        new_col_names = ["column%02d" % (c+1) for c in range(max(ncol_list))]
        common_header_keys['columns'] = dict(zip(new_col_names, [{}]*len(new_col_names)))
        
        with open(filename, 'w', newline='') as csvfile:
            csv_dialect = csv.excel()
            csv_dialect.delimiter = ";"
            writer = csv.writer(csvfile, dialect=csv_dialect)
            # First layer of keys
            row = []
            for key in common_header_keys:
                row += [key]+['']*np.max([0,len(common_header_keys[key])-1])
            writer.writerow(row)
            # Second layer of keys
            row = []
            for key in common_header_keys:
                if len(common_header_keys[key]) == 0:
                    row += [key]
                else:
                    row += [kkey for kkey in common_header_keys[key].keys()]
            writer.writerow(row)
            # Fill in values for each header
            for header in header_list:
                row = []
                for key in common_header_keys.keys():
                    if not key in header:
                        row += ['']*np.max([1,len(common_header_keys[key])])
                    elif len(common_header_keys[key].keys()) == 0:
                        row += [header[key]]
                    elif key == 'columns':
                        for keyi in range(len(common_header_keys[key])):
                            key_list = list(header['columns'].keys())
                            value_list = list(header['columns'].values())
                            if keyi < len(header['columns']):
                                row += [key_list[keyi]+": "+value_list[keyi]]
                            else:
                                row += ['']
                    else:
                        for kkey in common_header_keys[key].keys():
                            if kkey in header[key].keys():
                                row += [header[key][kkey]]
                            else:
                                row += ['']
                writer.writerow(row)

    def __dump_fasta(self, filename, header_list):
        sequences = []
        for header in header_list:
            file_id = header['filename'].rsplit(".",1)[0] + "_v" + str(header['version'])
            seq = Seq.Seq(header['protein']['sequence'])
            sequences.append(SeqRecord.SeqRecord(seq, id=file_id, description=""))
        with open(filename, "w", newline='') as file_handle:
            count = SeqIO.write(sequences, file_handle, "fasta")
        return(count)
            
class PrismData:
    """Base class of PRISM data types

    This class is thought as abstract and should not be instantiated. See documentation of derived classes.
    """
    
    def __init__(self, metadata, dataframe):
        """Constructor

        See documentation of derived classes
        """
        self.metadata = metadata
        self.dataframe = dataframe
        self.index_column_names = None

    def __repr__(self):
        """Representation - used e.g. for interactive printing"""
        return(self.header_to_string() + str(self.dataframe))

    def __str__(self):
        """Cast to string - used e.g. by print"""
        return("%s object of protein %s (uniprot %s) with dataframe shape %s" %
               (type(self), self.metadata['protein']['name'], self.metadata['protein']['uniprot'], str(self.dataframe.shape)))

    def copy(self):
        """Make a deep copy of self"""
        new_meta = copy.deepcopy(self.metadata)
        new_dataframe = self.dataframe.copy(deep=True)
        new_obj = PrismData(new_meta, new_dataframe)
        new_obj.__class__ = self.__class__
        return(new_obj)
        
    def seq_from_data(self):
        """Builds a sequence based on aa_ref and resi data columns

        Data is assumed to have index columns resi and aa_ref
        """

        if not 'aa_ref' in self.dataframe.columns:
            self.add_index_columns()
            
        n_res = 0
        seq = np.full(100,'X')
        for (aa_list,resi_list) in zip(self.dataframe['aa_ref'],self.dataframe['resi']):
            for (aa,resi) in zip(aa_list,resi_list):
                n_res = max(n_res,resi)
                if n_res > len(seq):
                    seq = np.concatenate([seq,np.full(max(100,n_res-len(seq)),'X')])
                if seq[resi-1] == 'X':
                    seq[resi-1] = aa
                elif seq[resi-1] != aa:
                    raise PrismFormatError("Data reference amino acid mismatch at position %d: % and %s" %
                                           (resi, seq[resi-1], aa))
        return("".join(seq[:n_res]))

    def check_header_and_data_sequence(self, verbose=0):
        """Check reference sequence from header against the data"""
        data_seq = self.seq_from_data()
        meta_seq = self.metadata['protein']['sequence']
        first_resn = int(self.metadata['protein']['first_residue_number']) if 'first_residue_number' in self.metadata['protein'].keys() else 1
        meta_seq = 'X'*(first_resn-1) + meta_seq # For checking purposes, matching indices makes things easier

        # Check sequence length
        if len(data_seq) < len(meta_seq):
            # Add X for unobserved positions
            data_seq = data_seq + 'X'*(len(meta_seq)-len(data_seq))
        elif (len(data_seq) != len(meta_seq)):
            raise PrismFormatError("Data has more residues than header sequence")
        
        # Check all substitution amino acids against header sequence
        match = [aa1==aa2 if aa1!='X' and aa2!='X' else True for (aa1,aa2) in zip(meta_seq,data_seq)]
        if not np.all(match):
            mismatch = np.where(np.logical_not(match))[0] + 1
            raise PrismFormatError("Mismatch between header sequence and data sequence at positions %s:\n%s\n%s" %
                                   (str(mismatch), meta_seq,data_seq))
        if (verbose > 0):
            print("Data and header amino acid sequence match with %d of %d positions not observed in data" % (data_seq.upper().count('X'),len(meta_seq)))
        return(data_seq)

    def header_to_string(self, skip_fields=[]):
        """Format meta-data into a string"""
        # Recursive function to write the next level of a dictionary
        def dump_subfields(dic, indent_level):
            return_str = ""
            dict_keys = list(dic.keys())
            # For the top level, order fields (dict is not ordered)
            if indent_level == 0:
                dict_keys_first = list()
                dict_keys_last = list()
                for k in ['version','protein']:
                    if k in dict_keys:
                        dict_keys.remove(k)
                        dict_keys_first.append(k)
                for k in ['variants','columns','filename']:
                    if k in dict_keys:
                        dict_keys.remove(k)
                        dict_keys_last.append(k)
                dict_keys = dict_keys_first + dict_keys + dict_keys_last
            # Make new indent and run recursively or dump values
            for key in dict_keys:
                if key in skip_fields:
                    continue
                if isinstance(dic[key], dict):
                    return_str += "# " + "    "*indent_level + key + ":\n"
                    return_str += dump_subfields(dic[key], indent_level+1)
                else:
                    # format floats values and use standard cast for others
                    s = "%.3f" % (dic[key]) if isinstance(dic[key],float) else str(str(dic[key])) 
                    return_str += "# " + "    "*indent_level + key + ": " + s + "\n"
            return(return_str)
        
        return_str = ""
        return_str += "# --------------------\n"
        return_str += dump_subfields(self.metadata, 0)                
        return_str += "# --------------------\n"
        return_str += "#\n"
        return(return_str)

    def check_column_names(self, verbose=0):
        # Check column names - case sensitive
        meta_colnames = list(self.metadata['columns'].keys()) + ['n_mut','aa_ref','resi','aa_var']
        data_colnames = self.dataframe.keys()[1:] # first column is data specific and not in header
        for cn in data_colnames:
            if not cn in meta_colnames:
                raise PrismFormatError("Could not find column name \'%s\' in header" % (cn))
        for cn in meta_colnames:
            if not cn in data_colnames:
                raise PrismFormatError("Could not find header column name \'%s\' in data" % (cn))
        if verbose > 0:
            print("Column names checks ok: %s" % (",".join(self.dataframe.keys())))
                

class VariantData(PrismData):
    """Container class for protein variant data

    Parameters
    ----------
    metadata : dict
        Metadata for data set
    dataframe : pandas.DataFrame
        Data

    Examples
    --------
    Get meta-data dictionary
    >>> meta_data = VariantData.metadata

    Get sequence from header
    >>> sequence = meta_data['protein']['sequence']

    Get pandas.DataFrame
    >>> data_frame = VariantData.dataframe

    Get all double mutants
    >>> data_frame[data_frame['n_mut']==2]

    Get all variants with substitutions into Ala
    >>> data_frame = VariantData.get_var_into_aa('A')

    Get all variants involving positions 30-45 (both included)
    >>> data_frame = VariantData.get_var_at_pos(range(30,45+1))

    Get a complete site-saturation (single mutant) library
    >>> data_frame = VariantData.saturate_variants()

    Merge three VariantData files and retrieve the union of variants (merge=outer)
    >>> merged_data = data1.merge([data2, data3], merge="outer")
    """
    # To be implmented:
    # Get all (multi-mutant) variants that involves D83V and L11P
    # >>> subst2var = VariantData.make_subst2var_map()
    # >>> data_frame = VariantData.dataframe[subst2var['D83V'] and subst2var['L11P']]

    def add_index_columns(self):
        var_split = [var.split(":") if var!="WT" else [] for var in self.dataframe['variant']]
        resi_ll = [[int(sub[1:-1]) for sub in var] if len(var)>0 else [] for var in var_split]
        aa_ref_ll = [[mut[0] for mut in var] if len(var)>0 else [] for var in var_split]
        aa_mut_ll = [[mut[-1] for mut in var] if len(var)>0 else [] for var in var_split]
        self.dataframe['n_mut'] = [len(var) for var in var_split]
        self.dataframe['aa_ref'] = aa_ref_ll
        self.dataframe['resi'] = resi_ll
        self.dataframe['aa_var'] = aa_mut_ll
        self.index_column_names = ['n_mut','aa_ref','resi','aa_var']

    def remove_index_columns(self):
        """Remove the indexing columns, i.e. ['n_mut','aa_ref','resi','aa_var']"""
        if 'aa_ref' in self.dataframe.columns:
            self.dataframe.drop(['n_mut','aa_ref','resi','aa_var'], axis=1, inplace=True)
        self.index_column_names = None

    def get_var_into_aa(self, aa, multimutant_mode='any'):
        """Get a data frame of variants with substitutions into a given amino acid

        Parameters
        ----------
        aa : character
            Amino acid in one-letter code
        multimutant_mode : string
            Mode of action for multi-mutants
            - any : At least one substitution per variant into aa required
            - all : All substitutions in variants must be into aa
            - exclude : Only return single mutants
        """
        # Check argument
        if not multimutant_mode in ['any','all','exclude']:
            raise ValueError("Function get_var_into_aa argument multimutant_mode must be \'any\', \'all\' or \'exclude\'")
        # Make a series of boleans to select rows
        if multimutant_mode=='any':
            bool_mask = [aa in aa_list for aa_list in self.dataframe['aa_var']]
        elif multimutant_mode=='all':
            bool_mask = [np.all(np.array(aa_list)==aa) if len(aa_list)>0 else False for aa_list in self.dataframe['aa_var']]
        elif multimutant_mode=='exclude':
            bool_mask = [aa_list==[aa] for aa_list in self.dataframe['aa_var']]
        else:
            assert False
        return(self.dataframe[bool_mask])
    
    def get_var_from_aa(self, aa, multimutant_mode='any'):
        # Exclude synonymous and nonsense?
        """Get a data frame of variants with substitutions from a given amino acid

        Parameters
        ----------
        aa : character
            Amino acid in one-letter code
        multimutant_mode : string
            Mode of action for multi-mutants
            - any : At least one substitution per variant from aa required
            - all : All substitutions in variants must be from aa
            - exclude : Only return single mutants
        """
        # Check argument
        if not multimutant_mode in ['any','all','exclude']:
            raise ValueError("Function get_var_from_aa argument multimutant_mode must be \'any\', \'all\' or \'exclude\'")
        # Make a series of boleans to select rows
        if multimutant_mode=='any':
            bool_mask = [aa in aa_list for aa_list in self.dataframe['aa_ref']]
        elif multimutant_mode=='all':
            bool_mask = [np.all(np.array(aa_list)==aa) if len(aa_list)>0 else False for aa_list in self.dataframe['aa_ref']]
        elif multimutant_mode=='exclude':
            bool_mask = [aa_list==[aa] for aa_list in self.dataframe['aa_ref']]
        else:
            assert False
        return(self.dataframe[bool_mask])
    
    def get_var_at_pos(self, target_resi, mode='any'):
        # Exclude synonymous and nonsense?
        """Get a data frame of variants with substitutions in one or more positions

        Parameters
        ----------
        target_resi : int or iterable of int
            Residue number as given in variants
        mode : string
            Mode of action for multi-mutants
            - any : At least one substitution in one of the given position
            - all : Substitutions at all given position
            - exact : Substitutions at all given position and no others
        """
        # Check argument mode
        if not mode in ['any','all','exact']:
            raise ValueError("Function get_var_from_aa argument mode must be \'any\', \'all\' or \'exact\'")
        # Check argument target_resi
        if isinstance(target_resi, int):
            target_resi = [target_resi]
        elif not isinstance(target_resi, list):
            target_resi = list(target_resi)
        target_resi = list(np.sort(target_resi))
        if mode == 'any':
            bool_mask = [np.any([ti in resi_list for ti in target_resi]) if len(resi_list)>0 else False for resi_list in self.dataframe['resi']]
        elif mode == 'all':
            bool_mask = [np.all([ti in resi_list for ti in target_resi]) if len(resi_list)>0 else False for resi_list in self.dataframe['resi']]            
        elif mode=='exact':
            bool_mask = [resi_list==target_resi for resi_list in self.dataframe['resi']]
        else:
            assert False
        return(self.dataframe[bool_mask])
            
    def get_subst(self, subst, single_only=True):
        """Get a data frame of variants that all contain a given substitution

        Parameters
        ----------
        subst : string
            Substitution, e.g. M1A
        single_only : bool
            Only return single mutants. If false, all multi-mutants with the given substitution will be returned
        """
        bool_mask = [subst in var.split(':') if var!='WT' else False for var in self.dataframe['variant']]
        # Add requirement on number of mutations if requested
        if single_only:
            bool_mask = [a and b for a,b in zip(bool_mask,self.dataframe['n_mut']==1)]
        # Return data frame slice
        return(self.dataframe[bool_mask])

    def order_variants(self, aa_order="CDEKRHNQAGSTVMLIFYWP=*~"):
        """Order variant according to number of substitutions, residue number, and a specified order of 
        amino acids. Non-specified amino acids are keps in original order. Inserts are not supported.
        """
        n_aa = len(aa_order)
        if len(set(aa_order)) != n_aa:
            raise ValueError("Argument aa_order to VariantData.order_variants must be unique characters")
        aa_dict = dict(zip(aa_order,range(n_aa)))
        aa_dict['WT'] = -1
        aa_list = [l[0] if len(l)>0 else 'WT' for l in self.dataframe['aa_var']]        
        self.dataframe['aa0'] = [aa_dict[aa] if aa in aa_dict.keys() else n_aa for aa in aa_list]
        self.dataframe['resi0'] = [l[0] if len(l)>0 else 0 for l in self.dataframe['resi']]
        self.dataframe.sort_values(by=['n_mut','resi0','aa0'], inplace=True)
        self.dataframe.index = range(self.dataframe.shape[0])
        self.dataframe.drop(columns=['aa0','resi0'], inplace=True)
        
    def saturate_variants(self, aa_order="CDEKRHNQAGSTVMLIFYWP", verbose=0):
        """Saturate variant list to make a full site-saturation library

        Order single mutants according to position and target amino acid and add NA 
        rows for missing single variants. Inserts are not supported.

        Parameters
        ----------
        aa_order : string
            Which amino acids to include and the order, default CDEKRHNQAGSTVMLIFYWP
        verbose : int, optional
            Level of output        
        """
        n_aa = len(aa_order)
        if len(set(aa_order)) != n_aa:
            raise ValueError("Argument aa_order to VariantData.saturate_variants must be unique characters")
        seq = self.metadata['protein']['sequence']
        first_resn = int(self.metadata['protein']['first_residue_number']) if 'first_residue_number' in self.metadata['protein'].keys() else 1
        var_list = [seq[resi]+str(resi+first_resn)+var if seq[resi] != var else seq[resi]+str(resi+first_resn)+"=" for resi in range(len(seq)) for var in aa_order]
        if verbose > 0:
            remove_set = set(self.dataframe['variant'])-set(var_list)
            if len(remove_set) > 0:
                if len(remove_set) > 500:
                    rm_list = list(remove_set)
                    print("Saturate variants will remove the following %d variants:\n%s\n... [skip %d variants] ...\n%s" %
                          (len(rm_list), ",".join([str(i) for i in rm_list[0:25]]), len(rm_list)-50, ",".join([str(i) for i in rm_list[-25:]])))
                else:
                    print("Saturate variants will remove the following %d variants:\n%s" % (len(remove_set), ",".join([str(i) for i in remove_set])))
        new_dataframe = pd.DataFrame({'variant':var_list})
        new_dataframe = new_dataframe.merge(self.dataframe, on='variant', how='left', copy=False)
        self.dataframe = new_dataframe
        self.add_index_columns()
        self.metadata['variants'] = self.calc_variants_metadata(verbose=verbose)

    def __missing_mask(self, verbose=0):
        """Private helper function for clean_na_variants and others
        
        Returns a bool array that marks rows with NA only. Should be fast."""
        # pd.DataFrame.iterrows() is not vectorized and very slow, DataFrame.apply is not too good either
        # return([np.all(row[1]) for row in self.dataframe.drop(columns=['variant','n_mut','aa_ref','resi','aa_var']).isnull().iterrows()])
        if 'aa_ref' in self.dataframe.columns:
            return(np.all(self.dataframe.drop(columns=['variant','n_mut','aa_ref','resi','aa_var']).isnull(), axis=1))
        else:
            return(np.all(self.dataframe.drop(columns='variant').isnull(), axis=1))

    def clean_na_variants(self, verbose=0):
        """Remove variants rows with all data missing"""
        missing_mask = self.__missing_mask()
        i = np.where(missing_mask)[0]
        if len(i) > 0:
            self.dataframe.drop(self.dataframe.index[i], axis=0, inplace=True)
            self.dataframe.index = range(self.dataframe.shape[0])
            if verbose > 0:
                print("Removed %d variant rows with all data missing" % (len(i)))

    def check_variants_metadata(self, no_overwrite=False, verbose=0):
        recalc_variants = self.calc_variants_metadata(verbose=verbose)
        if 'variants' in self.metadata:
            # Number of variants
            if 'number' in self.metadata['variants'].keys():
                if recalc_variants['number'] != self.metadata['variants']['number']:
                    s = "Header field \'variants:number\' is %d but %d variants (possible including WT) were found" % \
                        (self.metadata['variants']['number'], recalc_variants['number'])
                    if no_overwrite:
                        raise PrismFormatError(s)
                    elif verbose > 0:
                        print("WARNING: "+s)
            elif verbose>0 and not no_overwrite:
                print("Will add variants:number: %d to header" % (recalc_variants['number']))
                
            # Coverage
            if 'coverage' in self.metadata['variants'].keys():
                if np.abs(recalc_variants['coverage'] - self.metadata['variants']['coverage']) > 0.01:
                    s = "Header field \'variants:coverage\' is %.2f but variants calculated to cover %.2f" % \
                        (self.metadata['variants']['coverage'], recalc_variants['coverage'])
                    if no_overwrite:
                        raise PrismFormatError(s)
                    elif verbose > 0:
                        print("WARNING: "+s)
            elif verbose>0 and not no_overwrite:
                print("Will add variants:coverage: %.2f to header" % (recalc_variants['coverage']))
                    
            # Depth
            if 'depth' in self.metadata['variants'].keys():
                if np.abs(recalc_variants['depth'] - self.metadata['variants']['depth']) > 0.01:
                    s = "Header field \'variants:depth\' is %.2f but variants calculated to depth %.2f" % \
                        (self.metadata['variants']['depth'], recalc_variants['depth'])
                    if no_overwrite:
                        raise PrismFormatError(s)
                    elif verbose > 0:
                        print("WARNING: "+s)
            elif verbose>0 and not no_overwrite:
                print("Will add variants:depth: %.2f to header" % (recalc_variants['depth']))
         
            # Helper function
            def strip_all(s):
                return(s.lower().replace('.', '').replace('-','').replace('+','').replace('/','').replace(' ',''))

            # Variant width
            if 'width' in self.metadata['variants'].keys():
                if recalc_variants['width'] == "single mutants":
                    if not strip_all(self.metadata['variants']['width']) in ['singlemutants','singlemutant','singlemut','single','singleonly']:
                        s = "Header field \'variants:width\' is \'%s\' but only single mutants were found - use \'single mutants\'" % \
                            (self.metadata['variants']['width'])
                        if no_overwrite:
                            raise PrismFormatError(s)
                        elif verbose > 0:
                            print("WARNING: "+s)
                elif recalc_variants['width'] == "multi mutants":
                    if not strip_all(self.metadata['variants']['width']) in ['multimutants','multimutant','multimut','multi']:
                        # multimut = [var_list[i] for i in np.where(n_mut > 2)[0]]
                        s =  "Header field \'variants:width\' is \'%s\' but the following multi-mutant were found - use \'multi mutants\'" % \
                            (self.metadata['variants']['width'])
                        if no_overwrite:
                            raise PrismFormatError(s)
                        elif verbose > 0:
                            print("WARNING: "+s)
                elif recalc_variants['width'] == "single and double mutants":
                    if not strip_all(self.metadata['variants']['width']) in ['singleanddoublemutants','singleanddoublemutant','singleanddouble']:
                        s = "Header field \'variants:width\' is \'%s\' but variants are single and double mutants - use \'single and double mutants\'" % \
                            (self.metadata['variants']['width'])
                        if no_overwrite:
                            raise PrismFormatError(s)
                        elif verbose > 0:
                            print("WARNING: "+s)
                else:
                    raise ValueError("Recalculated variants width has an unknown value \'%s\'" % (recalc_variants['width']))
            elif verbose>0 and not no_overwrite:
                print("Will add variants:width: %s to header" % (recalc_variants['width']))
        else:
            if no_overwrite:
                raise PrismFormatError("No variants field in meta data")
                
        if not no_overwrite:
            self.metadata['variants'] = recalc_variants
        
    def calc_variants_metadata(self, verbose=0, return_missing_mask=False):
        ret = {}
        
        if verbose > 1:
            print("===============================================================")
            print("====    Calculate variants meta-data fields with timing    ====")
            start_time = time.process_time() 
            print("Calc number of variants (incl. missing_mask) ...", end=" ", flush=True)

        # Number of variants
        missing_mask = self.__missing_mask(verbose=verbose)
        ret['number'] = len(self.dataframe['variant']) - np.sum(missing_mask)

        if verbose > 1:
            elapsed_time = time.process_time() - start_time
            print("done. Elapsed time %.1f sec" % (elapsed_time))
            start_time = time.process_time() 
            print("Calc coverage and depth of variants ...", end=" ", flush=True)

        # Coverage and depth
        # Positions with substitutions
        non_missing_indices = np.where([not b for b in missing_mask])[0]
        flat_resi_list = [resi for resi_list in self.dataframe['resi'][non_missing_indices] for resi in resi_list]
        (resi_list,resi_counts) = np.unique(flat_resi_list, return_counts=True)
        
        # Coverage is positions with substitutions per residue in reference sequence
        ret['coverage'] = len(resi_list)/len(self.metadata['protein']['sequence'])
        # Depth is the average number of substitutions per position with substitutions
        ret['depth'] = np.mean(resi_counts)

        if verbose > 1:
            elapsed_time = time.process_time() - start_time
            print("done. Elapsed time %.1f sec" % (elapsed_time))
            start_time = time.process_time() 
            print("Calc width of variants ...", end=" ", flush=True)

        # Width
        n_mut = self.dataframe['n_mut']
        if np.all(n_mut <= 1):
            ret['width'] = "single mutants"
        elif np.max(n_mut) > 2:
            ret['width'] = "multi mutants"
        else:
            ret['width'] = "single and double mutants"

        if verbose > 1:
            elapsed_time = time.process_time() - start_time
            print("done. Elapsed time %.1f sec" % (elapsed_time))
            print("===============================================================")

        if return_missing_mask:
            return((ret,missing_mask))
        else:
            return(ret)
        
    def check(self, no_overwrite_metadata=False, verbose=0):
        """Perform all checks of data and meta data

        Checks formatting of all variants in the first column and sorts multimutants 
        and that the list is unique. Adds index columns.

        Consistancy checks of data and meta data includes column names, reference amino 
        acid sequence and all variants sub fields.

        Parameters
        ----------
        no_overwrite_metadata : bool
            Allow overwriting of the variants meta-data fields
        verbose : int, optional
            Level of output        
        """
        # Check and update variant column
        if verbose > 1:
            print("===============================================================")
            print("====     Check data section variant list with timing       ====")
            start_time = time.process_time() 
            print("Make var_list ...", end=" ", flush=True)

        var_list = [s.upper() for s in self.dataframe['variant']]

        # At least 2 characters for 'WT'
        if np.any([len(var) < 2 for var in var_list]):
            il = np.where([len(var) for var in var_list] < 2)[0]
            raise PrismFormatError("Bad variants: %s" % (str([var_list[i] for i in il])))

        if verbose > 1:
            elapsed_time = time.process_time() - start_time
            print("done. Elapsed time %.1f sec" % (elapsed_time))
            start_time = time.process_time() 
            print("Sort variants ...", end=" ", flush=True)

        # Sort multi-mutants according to residue number
        var_split_unordered = [var.split(":") if var!="WT" else [] for var in var_list]
        resi_ll_unordered = [[int(subst[1:-1]) for subst in var] if len(var)>0 else [] for var in var_split_unordered]

        # Make a new variant list with ordered multi-mutants
        var_split = []
        for ivar in range(len(var_split_unordered)):
            if len(set(resi_ll_unordered[ivar])) < len(resi_ll_unordered[ivar]):
                raise PrismFormatError("Variant %d %s has more substitutions on the same position" % (ivar+1, var_list[ivar]))
            var_split.append( [var_split_unordered[ivar][i] for i in np.array(resi_ll_unordered[ivar]).argsort()] )

        if verbose > 1:
            elapsed_time = time.process_time() - start_time
            print("done. Elapsed time %.1f sec" % (elapsed_time))
            start_time = time.process_time() 
            print("Split all substitutions ...", end=" ", flush=True)

        # Split substitutions in amino acids and residue number
        resi_ll = [[int(sub[1:-1]) for sub in var] if len(var)>0 else [] for var in var_split]
        aa_ref_ll = [[mut[0] for mut in var] if len(var)>0 else [] for var in var_split]
        aa_mut_ll = [[mut[-1] for mut in var] if len(var)>0 else [] for var in var_split]

        if verbose > 1:
            elapsed_time = time.process_time() - start_time
            print("done. Elapsed time %.1f sec" % (elapsed_time))
            start_time = time.process_time() 
            print("Make new_var_list ...", end=" ", flush=True)

        # Check all substitutions and make new variant list
        new_var_list = []
        for ivar in range(len(resi_ll)):
            subst_list = []
            for isub in range(len(resi_ll[ivar])):
                aa1 = aa_ref_ll[ivar][isub]
                resi = resi_ll[ivar][isub]
                aa2 = aa_mut_ll[ivar][isub]
                if not PrismParser.is_aa_one_nat(None, aa2, additional="*=~"):
                    raise PrismFormatError("Substitution %s%d%s is not to a standard amino acid" % (aa1,resi,aa2))
                if aa1 == aa2:
                    if verbose > 0:
                        print("Synonymous substitution %s%d%s renamed to %s%d=" % (aa1,resi,aa2,aa1,resi))
                    aa_mut_ll[ivar][isub] = '='
                subst_list.append("%s%d%s" % (aa_ref_ll[ivar][isub], resi_ll[ivar][isub], aa_mut_ll[ivar][isub]))
            new_var_list.append(":".join(subst_list) if len(subst_list)>0 else "WT")

        # Check if variants are unique
        (var_unique,var_counts) = np.unique(new_var_list, return_counts=True)
        if np.any(var_counts > 1):
            raise PrismFormatError("The following variants are listed multiple times:\n%s" % (var_unique[np.where(var_counts > 1)[0]]))

        # Update data frame
        self.dataframe['variant'] = new_var_list

        # Update index columns now that we have them
        self.dataframe['n_mut'] = [len(var) for var in var_split]
        self.dataframe['aa_ref'] = aa_ref_ll
        self.dataframe['resi'] = resi_ll
        self.dataframe['aa_var'] = aa_mut_ll
        self.index_column_names = ['n_mut','aa_ref','resi','aa_var']

        if verbose > 1:
            elapsed_time = time.process_time() - start_time
            print("done. Elapsed time %.1f sec" % (elapsed_time))
            print("===============================================================")

        self.check_column_names(verbose=verbose)
        self.check_header_and_data_sequence(verbose=verbose)
        self.check_variants_metadata(no_overwrite=no_overwrite_metadata, verbose=verbose)

    def to_new_reference(self, target_seq=None, first_resn=None, verbose=0,
                         min_identity=.8, min_coverage=.1, mismatch="remove", allow_inserts=True, allow_deletions=True):
        """Map variants from data file onto a new reference sequence
        
        This is the workhorse of data set merging

        Parameters
        ----------
        target_seq : string
            Amino acid sequence that variants are aligned to
        first_resn: integer >= 0
            If None, the header value is used if present, else one. If given the header value is replaced
        min_identity : float in [0.0, 1.0], optional
            Raise exception if the sequence identity is less than this (in covered region)
        min_coverage : float in [0.0, 1.0], optional
            Raise exception if the fraction of matched positions (i.e. not indel/gap) is less than this 
            This is a symmetric kind of coverage: (N_res-N_indel)/N_res where N_res is the lengths of the 
            longest sequence
        mismatch : string, optional
            How to handle mismatches in the alignment
            -convert: Change reference amino acid to match the new reference
            -remove: Remove variants that involve mismatch positions
            -not_allowed: Raise exception is a mismatch is incountered
        allow_inserts : bool, optional
            How to handle target sequence position missing in the sequence of the data 
            If true, the residue numbering will be ajusted accordingly
        allow_deletions : bool, optional
            How to handle data sequence positions missing in the target sequence
            If true, the residue numbering will be ajusted accordingly
        verbose : int, optional
            Level of output
        """

        # Determine new sequence offset and how much to shift residue numbering
        first_resn_meta = int(self.metadata['protein']['first_residue_number']) if 'first_residue_number' in self.metadata['protein'].keys() else 1
        if target_seq is None and first_resn is None:
            # Nothing to do
            return self
        elif first_resn is None:
            # Only alignment to target seq
            first_resn = first_resn_meta
            resi_shift_init = 0
        else:
            # Shift residue numbering (and possible align to target sequence)
            resi_shift_init = first_resn - first_resn_meta
            
        if first_resn < 0:
            raise ValueError("Argument first_resn to VariantData.to_new_reference must be a non-negative integer")

        seq = self.metadata['protein']['sequence']
        n_res_data = len(seq)
        
        # Which original data residue indices (0-based) to change
        resi_rm = []
        resi_shift = np.full(n_res_data,resi_shift_init)
        aa_change = {}
        if not target_seq is None:
            if not PrismParser.is_aa_one_nat(None, target_seq, "X"):
                raise ValueError("Argument target_seq to VariantData.to_new_reference must be a single-letter amino acid string (or None)")
            n_res_target = len(target_seq)
            # Align sequences
            # MatrixInfo supports a wildcard 'X' so use upper case. Open and extend penalty 11 and 1 is blast default
            # 2DO: Consider to write gap_function where X--XX is 2 openings and --XXX is only one, i.e. NT-gap is not an opening
            align = pairwise2.align.globalds(target_seq.upper(), seq.upper(), MatrixInfo.blosum62, -3, -1)
            # Check for alternative alignments
            if len(align) > 1 and verbose > 0:
                print("Found %d alignments" % (len(align)))
            if verbose > 0:
                # Dump verbose number of alignments
                for ia in range(min([verbose,len(align)])):
                    print("==== alignment %d of %d ====" % (ia+1,len(align)))
                    print(pairwise2.format_alignment(*align[ia]))
            # Check for equally good alignments, a good alignment has high score
            if len(align) > 1:
                align_scores = np.array([a[2] for a in align])
                if any(align_scores[0] - align_scores[1:] < 1e-3):
                    print("WARNING: There are alignments for %s with same or better scores (try verbose). Sores: %s" %
                          (self.metadata['protein']['name'],str(align_scores)))
            # Use first alignment
            align = align[0]
            assert len(align[0]) == len(align[1])
            n_align = len(align[0])
            # Number to convert 0-based alignment index to 0-based data (original) residue numbering
            ialign_shift = 0
            n_indel = n_match = n_mismatch = 0
            # Iterate over alignment, possible including gaps
            for ialign in range(len(align[0])):
                if align[0][ialign] == '-':
                    # Deletion: Position in data sequence not present in target sequence
                    if not allow_deletions:
                        raise PrismValueFail("Aborting alignment due to deletion at position %d which is not allowed" % (ialign+ialign_shift))
                    # Remove all substitutions at this position
                    resi_rm.append(ialign+ialign_shift)
                    # Downstream position are shifted to the left
                    resi_shift[(ialign+ialign_shift):] -= 1
                    n_indel += 1
                elif align[1][ialign] == '-':
                    # Insertion: Position in target sequence not present in data sequence
                    if not allow_inserts:
                        raise PrismValueFail("Aborting alignment due to insertion at position %d which is not allowed" % (ialign+ialign_shift))
                    resi_shift[(ialign+ialign_shift):] += 1
                    # Shift conversion number
                    ialign_shift -= 1
                    n_indel += 1
                elif align[0][ialign] != align[1][ialign]:
                    # Mismatch of 'X' always allowed since it is always mismatched
                    # Variants matched to X are always removed because position is ambiguous if adjacent to gap and because X is used for numbering offset
                    if align[0][ialign] == 'X' or align[1][ialign] == 'X':
                        resi_rm.append(ialign+ialign_shift)
                    elif mismatch == "not_allowed":
                        raise PrismValueFail("Aborting alignment due to mismatch at position %d which is not allowed" % (ialign+ialign_shift))
                    elif mismatch == "remove":
                        # Mark the substitutions at this position for removal
                        resi_rm.append(ialign+ialign_shift)
                    elif mismatch == "convert":
                        # Change reference amino acid
                        aa_change[ialign+ialign_shift] = align[0][ialign]
                    else:
                        raise ValueError("mismatch argument must be in [not_allowed,convert,remove]")
                    n_mismatch += 1
                else:
                    n_match += 1
            # Is alignment ok?
            n_res = max(n_res_target, n_res_data)
            coverage = (n_res-n_indel)/n_res
            if coverage < min_coverage:
                raise PrismValueFail("Sequence coverage %.3f (%d/%d) below threshold of %.3f" % (coverage, n_res-n_indel, n_res, min_coverage))
            identity = (n_res_target-n_mismatch)/n_res_target
            if identity < min_identity:
                raise PrismValueFail("Sequence identity %.3f below threshold of %.3f" % (identity, min_identity))
            # Update header if alignment is ok
            self.metadata['protein']['sequence'] = target_seq
            # Report actions
            if verbose > 0:
                print("Making new reference of %s with coverage %.3f and identity %.3f. Maximum shift in residue number is %d" %
                      (self.metadata['filename'], coverage, identity, max(abs(resi_shift))))
            if verbose > 0 and len(resi_rm) > 0:
                print("New reference sequence will remove variants at %d positions (mismatch setting: %s; indels: %d):\n%s" %
                      (len(resi_rm), mismatch, n_indel, ",".join([str(resi+first_resn_meta) for resi in resi_rm])))
            if verbose > 0 and len(aa_change) > 0:
                print("New reference sequence will change variants at %d positions (mismatch setting: %s; indels: %d):\n%s" %
                      (len(aa_change), mismatch, n_indel, ",".join([str(int(resi)+first_resn_meta) for resi in aa_change.keys()])))

        if not first_resn is None:
            self.metadata['protein']['first_residue_number'] = first_resn
        assert int(self.metadata['protein']['first_residue_number']) == first_resn
                
        # Split variant in amino acid letters and residue numbers
        aa_ref_ll = list(self.dataframe['aa_ref'])
        resi_ll = list(self.dataframe['resi'])
        aa_mut_ll = list(self.dataframe['aa_var'])

        # Update all variants
        new_var_list = []
        ivar_rm = []
        for ivar in range(len(resi_ll)):
            subst_list = []
            for isub in range(len(resi_ll[ivar])):
                # 0-based index of original sequence without offset
                resi_orig = resi_ll[ivar][isub] - first_resn_meta
                # Is substitution marked for removal?
                if resi_orig in resi_rm:
                    ivar_rm.append(ivar)
                    break
                # Update residue number
                resi_ll[ivar][isub] += resi_shift[resi_orig]
                # Update reference amino acid
                if resi_orig in aa_change.keys():
                    if aa_change[resi_orig] == aa_mut_ll[ivar][isub]:
                        if verbose > 2:
                            print("Data contains a substitution, %s%d%s, into the new reference which will be removed" %
                                  (aa_ref_ll[ivar][isub],resi_ll[ivar][isub],aa_mut_ll[ivar][isub]))
                        ivar_rm.append(ivar)
                        break
                    else:
                        if verbose > 2:
                            print("Convert variant due to mismatch: %s%d%s -> %s%d%s" %
                                  (aa_ref_ll[ivar][isub], resi_ll[ivar][isub]-resi_shift[resi_orig], aa_mut_ll[ivar][isub],
                                   aa_change[resi_orig], resi_ll[ivar][isub], aa_mut_ll[ivar][isub]))
                        aa_ref_ll[ivar][isub] = aa_change[resi_orig]
                subst_list.append("%s%d%s" % (aa_ref_ll[ivar][isub], resi_ll[ivar][isub], aa_mut_ll[ivar][isub]))
            if len(subst_list) > 0:
                new_var_list.append(":".join(subst_list))
            else:
                new_var_list.append("WT")
                
        if verbose > 0 and len(ivar_rm) > 0:
            rm_list = [str(self.dataframe['variant'].iat[i]) for i in ivar_rm]
            if len(rm_list) > 500:
                print("New reference sequence will remove %d variants (mismatch setting: %s; indels: %d):\n%s\n... [skip %d variants] ...\n%s" %
                          (len(rm_list), mismatch, n_indel, ",".join(rm_list[0:25]), len(rm_list)-50, ",".join(rm_list[-25:])))
            else:
                print("New reference sequence will remove %d variants (mismatch setting: %s; indels: %d):\n%s" %
                          (len(rm_list), mismatch, n_indel, ",".join(rm_list)))

        self.dataframe['variant'] = new_var_list
        self.dataframe['aa_ref'] = aa_ref_ll
        self.dataframe['resi'] = resi_ll
        self.dataframe['aa_var'] = aa_mut_ll
        # Remove rows containing variants with deleted substitutions
        self.dataframe = self.dataframe.drop(self.dataframe.index[ivar_rm], axis=0)
        self.dataframe.index = range(self.dataframe.shape[0])

    def merge(self, data_list, target_seq=None, first_resn=None, merge='outer', **kwargs):
        """Merge a list of data files into a single DataFrame and dump it in a file

        Parameters
        ----------
        data_list : list of VariantData
            The list of data sets to merge
        target_seq : string
            Amino acid sequence that all data sets are aligned to. If None, 
            the sequence of the first data set will be used.
        first_resn : Integer >= 0
            Residue number of first amino acid in target_seq, default one. If 
            target_seq is None, the header value will be used if present or 
            overwritten if an argument is given here.
        merge : string
            How to merge:
            - outer: output the union of variants, default
            - left: Only keep veriants that are present in the self data set
            - inner: output the intersect of variants
        **kwargs : keyword arguments
            Passed to to_new_reference function
        """
        if not merge in ['left', 'outer', 'inner']:
            raise ValueError("Allowed merge arguments are left, outer or inner")

        merged_data = self.copy()

        if isinstance(data_list, VariantData):
            data_list = [data_list]
        elif len(data_list) < 1:
            raise ValueError("List of data set to merge is empty")
        
        # Target sequence
        if target_seq is None:
            # Make from meta data, variant residue numbers will match the index of this
            target_seq = self.metadata['protein']['sequence']
            if not first_resn is None:
                raise ValueError("merge argument first_resn can only be set if target_seq != None\n" + \
                                 "Use VariantData.to_new_reference to only shift residue numbering")
            first_resn = int(self.metadata['protein']['first_residue_number']) if 'first_residue_number' in self.metadata['protein'].keys() else 1
        else:
            # Align self data set to target
            if first_resn is None:
                first_resn = 1
            else:
                if first_resn < 0:
                    raise ValueError("Argument first_resn must be non-negative")
            try:
                merged_data.to_new_reference(target_seq, first_resn, **kwargs)
            except PrismValueFail as msg:
                print("Could not align data to target sequence: %s" % (msg))
                return(None)
            
        # Remove the index columns ['n_mut','aa_ref','resi','aa_mut'] before merging
        merged_data.remove_index_columns()
        # Add '_00' to the names of all but variant column
        merged_data.dataframe.columns = [merged_data.dataframe.columns[0]] + [s+"_00" for s in merged_data.dataframe.columns[1:]]

        # Init new meta data dictionary
        merged_metadata = {'version':0, 'protein':{}, 'merged':{}, 'variants':{}, 'columns':{}}
        merged_metadata['protein'] = copy.deepcopy(self.metadata['protein'])
        
        merged_metadata['merged']["file_00"] = self.metadata['filename']+" (version "+str(self.metadata['version'])+")"
        for key in self.metadata['columns'].keys():
            merged_metadata['columns'][key+"_00"] = self.metadata['columns'][key]

        # Align and merge other data sets
        c = 0
        for data in data_list:
            data_copy = data.copy()
            # 2DO: Consider to check uniprot (others are human readable and may have typoes) of data to be merged
            try:
                data_copy.to_new_reference(target_seq, first_resn, **kwargs)
            except PrismValueFail as msg:
                print("Skip merge of file %s: %s" % (data_copy.metadata['filename'],msg))
                continue
            data_copy.remove_index_columns()
            suffix = "_%02d" % (c+1)
            data_copy.dataframe.columns = [data_copy.dataframe.columns[0]] + [s+suffix for s in data_copy.dataframe.columns[1:]]
            merged_data.dataframe = merged_data.dataframe.merge(data_copy.dataframe, on='variant', how=merge, suffixes=(False,False), copy=False)
            if merged_metadata['protein']['uniprot'] != data.metadata['protein']['uniprot']:
                print("WARNING: Uniprot mismatch on merge: Merging file %s uniprot %s with file %s of different uniprot %s" % (
                    merged_metadata['merged']["file_00"], merged_metadata['protein']['uniprot'],
                    data.metadata['filename'], data.metadata['protein']['uniprot']))
            merged_metadata['merged']["file_%02d" % (c+1)] = data_copy.metadata['filename']+" (version "+str(data_copy.metadata['version'])+")"
            for key in data_copy.metadata['columns'].keys():
                merged_metadata['columns'][key+"_%02d" % (c+1)] = data_copy.metadata['columns'][key]
            c += 1
                
        # Add the index columns ['n_mut','aa_ref','resi','aa_mut'] again
        merged_data.add_index_columns()
        
        merged_metadata['protein']['sequence'] = target_seq
        if first_resn != 1:
            merged_metadata['protein']['first_residue_number'] = first_resn
        merged_data.metadata = merged_metadata
        merged_data.metadata['variants'] = merged_data.calc_variants_metadata()
        return(merged_data)

        
# class ResidueData(PrismData):
#     """Container class for per-residue data
#     """
#
#     def _init(self):
#         self.index_column_names = ['aa_ref']
#
#     def add_index_columns(self):
#         pass
#
#     def remove_index_columns(self):
#         """Remove the indexing columns, i.e. 'aa_ref'"""
#         self.dataframe.drop(['aa_ref'], axis=1, inplace=True)
#
#     def check(self, verbose=0):
#         # check ResidueData specific stuff
#         check_status = None
#         try:
#             self.check_column_names(verbose=verbose)
#             self.check_header_and_data_sequence(verbose=verbose)
#             self.check_variants_metadata(no_overwrite=no_overwrite_metadata, verbose=verbose)
#         except PrismFormatError as msg:
#             raise PrismFormatError("VariantData.check failed. First fail message:\n%s" % (msg))
#         else:
#             check_status = True
#
#         return(check_status)
#
#     def to_new_reference(self, target_seq=None, first_resn=None, verbose=0):
#         pass


if __name__ == "__main__":
    # Parse commandline arguments
    import argparse,sys
    arg_parser = argparse.ArgumentParser(description="PRISM data file processing and alignment")
    # positional arguments
    arg_parser.add_argument("files", nargs="+", metavar="DATAFILE",
                        help="Data files to process")
    # keyword arguments
    arg_parser.add_argument("--dump_header_csv", metavar="FILE",
                      help="Dump a .csv file containing an overview of file headers")
    arg_parser.add_argument("--dump_fasta", metavar="FILE", 
                      help="Dump all sequences in a fasta file")
    arg_parser.add_argument("--parse", action='store_true',
                      help="Read data file(s), run all checks and write new data file(s)")
    arg_parser.add_argument("--merge", metavar="FILE", 
                      help="Merge data files concerning the (approximately) same target protein")
    # A run of this main can only have a single target sequence so to parse more files with different sequence, run multiple times
    arg_parser.add_argument("--target_seq", metavar="FILE_OR_SEQUENCE",
                      help="Align all output to this sequence, default is the sequence given in the first file")
    arg_parser.add_argument("--saturate", action='store_true',
                      help="In the output data, add rows to cover all 20 substitutions at all considered positions")
    arg_parser.add_argument("--clean_na", action='store_true',
                      help="In the output data, add rows to cover all 20 substitutions at all considered positions")
    arg_parser.add_argument("--first_res_num", metavar="INT", 
                      help="In the output data, assign this number to the first given amino acid of the sequence")
    arg_parser.add_argument("-v","--verbose", action="count", default=0,
                      help="Level of output, default zero is no output")
    args = arg_parser.parse_args()
    if args.saturate and args.clean_na:
        raise ValueError("Arguments saturate and clean_na cannot be used together")
    
    # Read data files
    data_parser = PrismParser()
    data_list = []
    for filename in args.files:
        try:
            data = None
            if args.merge or args.parse:
                if args.verbose > 0:
                    print("\nRead data file %s " % (filename),flush=True)
                data = data_parser.read(filename, verbose=args.verbose)
            elif args.dump_header_csv or args.dump_fasta:
                # Only read header which is much faster for large files
                if args.verbose > 0:
                    print("\nRead header of %s " % (filename),flush=True) 
                header = data_parser.read_header(filename, verbose=args.verbose)
                data = PrismData(header, None)
            else:
                break
        except PrismFormatError as error:
            print("ERROR reading file %s:\nBad file format: %s\n" % (filename,error))
        except yaml.parser.ParserError as error:
            print("ERROR reading file %s:\nBad YAML syntax (ParserError) in header:\n%s\n" % (filename,error))
        except yaml.scanner.ScannerError as error:
            print("ERROR reading file %s:\nBad YAML syntax (ScannerError, perhaps a missing whitespace after colon?) in header:\n%s\n" % (filename,error))
        except FileNotFoundError:
            print("ERROR reading file %s:\nCannot find file\n" % (filename))
        else:
            # Only if all is good
            data_list.append(data)

    # Check input
    if len(data_list) < 1:
        print("\nERROR: Missing data files or arguments")
        arg_parser.print_usage()
        sys.exit(2)
    if args.verbose > 0:
        print("Done loading %d file(s)\n" % (len(data_list)))        

    # Dump overview csv if requested
    if args.dump_header_csv:
        if args.target_seq or args.first_res_num or args.saturate or args.clean_na:
            print("WARNING: Arguments target_seq, first_res_num, saturate and clean_na have no effect on dump_header_csv")
        data_parser.write_meta_file(data_list, args.dump_header_csv, output_type="csv")
        print("Dumped header overveiw of %d files in %s" % (len(data_list),args.dump_header_csv))

    # Dump amino acid sequences in fasta format if requested
    if args.dump_fasta:
        if args.target_seq or args.first_res_num or args.saturate or args.clean_na:
            print("WARNING: Arguments target_seq, first_res_num, saturate and clean_na have no effect on dump_fasta")
        data_parser.write_meta_file(data_list, args.dump_fasta, output_type="fasta")
        print("Dumped %d sequences in %s" % (len(data_list),args.dump_fasta))

    # Set target sequence and offset for all data files
    target_seq = None
    first_resn = 1
    if args.target_seq:
        if args.target_seq[-6:].lower() == ".fasta" or args.target_seq[-4:].lower() == ".fas":
            record = None
            for r in SeqIO.parse(args.target_seq, "fasta"):
                if not record is None:
                    # if args.verbose > 0:
                    print("WARNING: Only using the first sequence record in %s" % (args.target_seq))
                    break
                record = r
            target_seq = str(record.seq)
            print("Target sequence \'%s\' of %d residues from file %s" % (record.description, len(target_seq), args.target_seq))
        else:
            target_seq = args.target_seq
        if args.first_res_num:
            first_resn = int(args.first_res_num)
    else:
        target_seq = (data_list[0]).metadata['protein']['sequence']
        if args.first_res_num:
            # Overwrite offset given in header of first file
            first_resn = int(args.first_res_num)
        elif 'first_residue_number' in (data_list[0]).metadata['protein']:
            # Only use offset of first file with that sequence (and not the commandline input)
            first_resn = int((data_list[0]).metadata['protein']['first_residue_number'])

    if not data_parser.is_aa_one_nat(target_seq, "X"):
        print("ERROR: Target sequence seems not to be one-letter amino acids or \'X\':")
        print("Target sequence: %s" % (target_seq))
        sys.exit(2)    

    # Parse & dump all files if requested
    if args.parse:
        filename_suffix = "_parsed"
        for data in data_list:
            if args.target_seq or first_resn:
                data.to_new_reference(target_seq, first_resn, verbose=args.verbose,
                                      min_identity=.9, min_coverage=.1, mismatch="remove", allow_inserts=True, allow_deletions=True)
            if args.saturate:
                data.saturate_variants()
            elif args.clean_na:
                data.clean_na_variants()
            filename_parts = data.metadata['filename'].split("/")[-1].rsplit(".",1)
            outfilename = filename_parts[0] + filename_suffix + "." + filename_parts[1]
            data_parser.write(outfilename, data)
        print("Parsed and wrote %d data file(s) with output filename suffix %s" % (len(data_list), filename_suffix))
    
    # Dump a merged data file if requested
    if args.merge:
        merged_data = data_list[0].merge(data_list[1:], target_seq, first_resn, merge='outer', verbose=args.verbose,
                                         min_identity=.9, min_coverage=.1, mismatch="remove", allow_inserts=True, allow_deletions=True)
        if merged_data is None:
            print("ERROR: Could not merge data")
            sys.exit(2)    
        
        if args.saturate:
            merged_data.saturate_variants()
        elif args.clean_na:
            merged_data.clean_na_variants()
        data_parser.write(args.merge, merged_data)
        print("Dumped merge of %d data files in %s " % (len(merged_data.metadata['merged']), args.merge))
