#!/bin/bash

if [[ $1 = '--Help' ]]; then

   bin/admevolve --Help
   exit

fi

build.sh admevolve.gpr || exit

rm -rf results
mkdir -p results
touch results/history.txt

bin/admevolve $* \
   --Courant 0.25 \
   --Tfinal 11.0 \
   --PrintCycle 10 \
   --PrintTimeStep 11.0 \
   --MaxTimeSteps 40000 \
   --NumTasks 8 \
   --OutputDir results \
   --UseRendezvous \
   --DataDir data | tee adm-history.txt

# bin/admevolve $* \
#    --Courant 0.25 \
#    --Tfinal 11.0 \
#    --PrintCycle 10 \
#    --PrintTimeStep 11.0 \
#    --MaxTimeSteps 40000 \
#    --NumTasks 8 \
#    --OutputDir results \
#    --UseProtObject \
#    --DataDir data | tee adm-history.txt

# bin/admevolve $* \
#    --Courant 0.25 \
#    --Tfinal 11.0 \
#    --PrintCycle 10 \
#    --PrintTimeStep 11.0 \
#    --MaxTimeSteps 40000 \
#    --NumTasks 8 \
#    --OutputDir results \
#    --UseTransientTasks \
#    --DataDir data | tee adm-history.txt

# bin/admevolve $* \
#    --Courant 0.25 \
#    --Tfinal 11.0 \
#    --PrintCycle 10 \
#    --PrintTimeStep 11.0 \
#    --MaxTimeSteps 40000 \
#    --NumTasks 8 \
#    --OutputDir results \
#    --UseSyncBarriers \
#    --DataDir data | tee adm-history.txt
