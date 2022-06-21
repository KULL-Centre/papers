Data and code for "Prediction of quality-control degradation signals in yeast proteins" by Johansson et al

The degron prediction program QCDpred is available as a stand-alone python program and a web server using Googles CoLab Service 

QCDpred.py
----------

The python script only depends on the standard modules os, sys, argparse and numpy. It has been tested with python 3.6 and 3.9
but will likely work with any python version >3.0. The model is fast and should be able to process a proteome in less than a minute.

Example of terminal usage
$ python3 QCDpred.py CALLQSRLLLSAPRRAAATARY
          seq00000  CALLQSRLLLSAPRRAA  0.87015  L     9
          seq00000  ALLQSRLLLSAPRRAAA  0.83666  L    10
          seq00000  LLQSRLLLSAPRRAAAT  0.80842  S    11
          seq00000  LQSRLLLSAPRRAAATA  0.73379  A    12
          seq00000  QSRLLLSAPRRAAATAR  0.61286  P    13
          seq00000  SRLLLSAPRRAAATARY  0.74813  R    14

Input either a sequence or FASTA file with one or more amino acid sequences. All sequences should be 17 amino acids or longer in whith case the script will tile the sequence into fragments of 17 amino acids. The scripts outputs one line per tile giveing the sequence identifier, the tile sequence, a degron probability, and the central amino acid identity and number. 

Webserver
---------

The file QCDpred.ipynb is made to run as a webservice at Google colab available at
https://colab.research.google.com/github/KULL-Centre/papers/blob/main/2022/degron-predict-Johansson-et-al/QCDpred.ipynb

This requires a login for Googles services, e.g. a gmail, and will enable you to store a copy of the program in your Google drive
