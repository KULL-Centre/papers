import pandas as pd
import os

class record(): 
    """A class to hold each sequencing item. Each item has a name,
    a total amount of mutations and a list of each mutation, that 
    has the format (G192T).
    """
    def __init__(self,name,total,mutList):
        self.name = name
        self.total = total
        self.mutList = mutList

def parse_mut_file(f):
    """This function parses a mut file that contains information 
    on all mutations different sequences in the following format:
    
    >name
    total 2
    2
    T 27 A
    G 192 T
    
    The function outputs a list of record objects for each entry
    in the given mutfile. A list is used to preserve the same 
    sequence of entries.
    
    (str) -> (list)
    """
    
    #Open the file in the path passed as an argument
    file = open(f,'r')
    
    #Container for the output
    records = []
    
    #The first line will be the name of the first entry
    name = next(file)[1:].strip()
    
    name = None
    total = None
    mutList = []    
    
    #Loop over the rest of the file
    for line in file:
        line = line.strip()
        
        #Each time another name is reached, the entry is added to
        #the output bin
        if line[0] == '>':
            #Save the previous record
            records.append(record(name,total,mutList))
            
            #Reset the parameters that were just saved
            name = None
            total = None
            mutList = []
            
            #Start creating a new record
            name = line[1:]
        
        #Check the total amount of mutations
        elif 'total ' in line:
            total = line[6:]
            
        else:
            #Check for the line with just the total number
            try:
                int(line)
            #The only other type of line in the file is a mutation
            except:
                old = line[0]
                new = line[-1]
                
                #Add the mutation to the mutList and change the format
                #from (A 12 T) to (A12T)
                mutList.append(''.join([old,
                                     line[1:-1].strip(),
                                     new,
                                     ' ']))
    
    #Save the last entry            
    records.append(record(name,total,mutList))
    
    return records
    
def bin_reads(l):
    """The function counts the number of entries in a list that
    have the same mutations. This information is the returned as
    a dictionary where the keys are all the sequences of mutations
    and the values are the number of entries with that given
    sequence of mutations. A dictionary is used because of
    improved performance when checking for a given key.
    
    (list) -> (dict)
    """
    
    #Dictionary to hold the results
    out = {}
    
    #Loop over all the entries in the list
    for read in l:
        #The mutations are combined to a single string
        muts = "".join(read.mutList).strip()
        
        #Check if this is in the dictionary and add 1 
        #to the count. Otherwise create a dictionary 
        #entry for this sequence of mutations
        if muts in out:
            out[muts] += 1
        else:
            out[muts] = 1

    return out
        
if __name__ == "__main__":
    #Loop over the files in the current directory
    for file_name in os.listdir():
        if '_AlignmentOutput_DNA.txt' in file_name:
            
            #Parse the mut textfile
            sample = parse_mut_file(file_name)
            
            #Count all the reads with same sequence
            sample_c = bin_reads(sample)
            
            #Create a pandas dataframe for easy saving
            out = pd.DataFrame(sample_c.items())
            
            #Save the dataframe as a csv file
            out.to_csv(''.join([file_name[:-4]
                                ,'_count.csv']))
    
