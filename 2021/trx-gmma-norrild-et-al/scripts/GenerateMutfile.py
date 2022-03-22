from Bio import SeqIO
from Bio.Seq import Seq
import os

#The template for the full sequenced region
template = Seq(''.join(['CCGTAAAGATAACGACGAAGAGGCCAAGAAGGTTGAGTATATT',
                        'GTGCGCGAACTGGCGCAGGAATTTGACGGTCTGATCATGGTTT',
                        'TCGAGCTGGACACGAACAAGGCACCGGAGATCGCGAAAAAGTA',
                        'CAATATCACCACCACCCCGACTGTCGCATTTTTCAAAAATGGC',
                        'GAGGTCAAGAGCGTTCTGATTGGCGCGATTCCAAAAGACCAGC',
                        'TGCGTGATGAAATCCTGAAATATCTGGGTCACCATCATCACCA',
                        'TCACGGTACCAAACCGTACCAACGTCAGTTCATCGA']))
  
#Slicing the template to the mutated region
template = template[71:217]

#Loop over all the fasta files in the current folder with
#with the merged sequences from each sample
for file_name in os.listdir():
    if '.fasta' in file_name:
        #Open the output file
        myfile = open(''.join([file_name[:-6], #This is the sample name in my format
                               '_AlignmentOutput_DNA.txt']),'w')
        
        #For checking progress of the code
        count = 0         
            
        #Looping over all fasta entries in the specified file
        for record in SeqIO.parse(file_name,'fasta'):
            count += 1
            
            #Only the mutated region of oligo 4 & 5 is to be used
            read = record.seq[71:217]
            
            # A list to hold the mutations in the format of strings "L 11 P"
            substitutions = []
            
            #Check for positions that differ and append to the substitutions list 
            for i in range(len(read)):
                if read[i] != template[i]:
                    substitutions.append(''.join([str(template[i]),
                                                  ' ',
                                                  str(i+71+1),#71 is the first nucleotide in the mutated region. When numbering from 1 instead of 0, it will be 72   
                                                  ' ',
                                                  str(read[i])]))                 
            
            #Write the entry to the results file
            myfile.write(''.join(['>',str(record.name),'\n']))
            myfile.write(''.join(['total ',str(len(substitutions)),'\n']))
            myfile.write(''.join([str(len(substitutions)),'\n']))
            for sub in substitutions:
                myfile.write(''.join([sub,'\n']))
            
            #Print counts for process 
            if count % 5000 == 0:
                print(count)
        
        myfile.close()