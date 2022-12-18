#!/bin/sh


## MAKE THE REQUIRED BINARIES
make nuke && make


## GENERATE THE REQUIRED DATA SETS
R CMD BATCH --no-save generate_data.R /dev/stdout


## RUN THE EXAMPLES (THIS WILL TAKE A WHILE!)
mkdir -p results/{gsgp,baggp,2segp,2segp-noboot,sspbe}

parallel --eta --bar --progress \
	 ./dist/gsgp ./problems/seeds/{1} {2} D ./problems/data/{1} ./problems/splits/{1} {2} 30 {3} 6 200 250 '>' ./results/gsgp/{1}-{2}-{3} \
	 ::: mbishop-01 mbishop-05 mbishop-10 ::: `seq 30` ::: 0.001 0.002 0.004 0.008 0.016 0.031 0.062 0.125 0.250 0.500 1.000

problem=friedman1
parallel --eta --bar --progress \
	 ./dist/gsgp ./problems/seeds/$problem {1} D ./problems/data/$problem ./problems/splits/$problem {1} 30 {2} {3} 200 250 250 n '>' ./results/gsgp/$problem-{1}-{2}-{3} \
	 ::: `seq 30` ::: 0.01 0.1 0.25 0.5 0.75 1.00 ::: `seq 20`

problem=friedman1
parallel --eta --bar --progress \
	 ./dist/gsgp ./problems/seeds/$problem {1} D ./problems/data/$problem ./problems/splits/$problem {1} 30 {2} 6 200 250 1 '>' ./results/gsgp/$problem-{1}-{2}-250gen \
	 ::: `seq 30` ::: 0.01 1.00

problem=friedman1
parallel --eta --bar --progress \
	 ./dist/baggp ./problems/seeds/$problem {1} ./problems/data/$problem ./problems/splits/$problem {1} 30 50 500 100 '>' ./results/baggp/$problem-{1} \
	 ::: `seq 30`

problem=friedman1
parallel --eta --bar --progress \
	 ./dist/2segp ./problems/seeds/$problem {1} D ./problems/data/$problem ./problems/splits/$problem {1} 30 {2} 500 100 '>' ./results/2segp/$problem-{1}-{2} \
	 ::: `seq 30` ::: `seq 4` `seq 5 5 50`

problem=friedman1
parallel --eta --bar --progress \
	 ./dist/2segp ./problems/seeds/$problem {1} D ./problems/data/$problem ./problems/splits/$problem {1} 30 -{2} 500 100 '>' ./results/2segp-noboot/$problem-{1}-{2} \
	 ::: `seq 30` ::: `seq 4` `seq 5 5 50`

problem=friedman1
parallel --eta --bar --progress \
	 ./dist/sspbe ./problems/seeds/$problem {1} D ./problems/data/$problem ./problems/splits/$problem {1} 30 {2} 192 250 '>' ./results/sspbe/$problem-{1}-{2} \
	 ::: `seq 30` ::: 1 2 3 5 8 13 21 34 55 89 144 196


## PERFORM THE DECOMPOSITION ON THE RESULTS (takes about 30 min on my machine)
R CMD BATCH --no-save decompose_results.R /dev/stdout


## CREATE THE FIGURES IN THE PAPER
R CMD BATCH --no-save generate_figures.R /dev/stdout
