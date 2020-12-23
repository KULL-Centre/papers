#!/bin/bash

for t in 50 150 250 400 1000
do
	nohup sh Flow_Backmapping.sh ${t} &
done
