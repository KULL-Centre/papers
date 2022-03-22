from Bio import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
import os 

#Reading the fastq files from the same directory
for file in os.listdir():
    if file.find('R1') != -1:
        file1 = file
    if file.find('R2') != -1:
        file2 = file

#List containing Seq.Record files to be written to the output file
with open('test.fasta','a') as myFile:
    
    secondFile = SeqIO.parse(file2,'fastq')
    
    count = 0
    succes = 0
    for read1 in SeqIO.parse(file1,'fastq'):
        count += 1
        read2 = next(secondFile)
        
        #Identifying the amplicon primer ends in both reads
        #There is one position that consistently is a N instead
        #of the correct G. I think that this must be a bad cycle
        #or some quirk of the sequencing.
        if ('ccgtaaagataacgacgaagaggc'.upper() in read1.seq and 
            ('tcgatgaactgacgttggtacgg'.upper() in read2.seq or 
             'tcgatgaantgacgttggtacgg'.upper() in read2.seq)):
            fw = read1.seq
            rev = read2.seq
            
        elif ('ccgtaaagataacgacgaagaggc'.upper() in read2.seq and 
              ('tcgatgaactgacgttggtacgg'.upper() in read1.seq or 
               'tcgatgaantgacgttggtacgg'.upper() in read1.seq)):
            fw = read2.seq
            rev = read1.seq
            
        else:
            fw = ''
            rev = ''
        
        #Check if the reads are identical in the mutated area
        if rev != '' and fw != '':
            if fw[61:227] == rev[67:233].reverse_complement():
                #Join the sequences outside of the mutated area
                #The parts in parenthesis are omitted 
                #  
                # --------|------------------>|(--->)
                #         |                   |
                #   (<---)|<------------------|-------------
                #
                cand = SeqRecord(str(fw[:227] + rev[:67].reverse_complement()),
                                 id = read1.id,
                                 name = read1.name,
                                 description=read1.description)
                
                #Correcting the quirk of the dataset
                dummy = MutableSeq(cand.seq)
                dummy[-9] = 'G'
                cand.seq = dummy.toseq()
                
                #Checking for the correct length
                if len(cand.seq) == 294:
                    
                    SeqIO.write(cand,myFile,'fasta')
                    succes += 1
                    
        #Printing the progress of the code               
        if count % 25000 == 0:
            print(count)
            print(succes)
    
    #Finalizing the script and saving the merged fasta files
    print('The total counts is: ' + str(count))
    print('The total succes is: ' + str(succes))


    
    