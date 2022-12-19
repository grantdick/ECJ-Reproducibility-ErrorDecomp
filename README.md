# Error Decomposition for Improved Reproducibility in Evolutionary Computation

Code required for paper submitted to [special issue of ECJ on reproducibility in evolutionary computation](https://direct.mit.edu/DocumentLibrary/CallForPapers/ecj-si-rep.pdf). Provides a demonstration of how error decomposition can lead to better understanding of algorithm behaviour (and therefore better reproduction and interpretation of results).

## Running the Experiments

The project should run under Mac OS or Linux (incl. WSL) without problem. GCC or clang is required for compilation, [GNU parallel](https://www.gnu.org/software/parallel/) is required to run the experiments, and R is required for data preparation and analysis. Within R, we use the [tidyverse](https://www.tidyverse.org/) and [ggthemes](https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/) libraries. Otherwise, the C code is standalone and does not require any additional libraries.  

A complete cycle of the project (compilation, running, and producing figures) can be done by running the script `perform_runs.sh`.

The C code is not particularly memory intensive, so multiple instances can be run simultaneously (as is done in the shell script using parallel) without too many worries. However, even using parallel processing, running the experiments will take some time (on a Threadripper 3960X, running the code takes about 15 hours). Performing the error decomposition on the results takes about 30 minutes on a single core (R seems to use a lot of RAM for this, so 32GB is recommended), and final visualisation of the results takes a few seconds.

## Questions/Help
Please [email me](mailto:grant.dick@otago.ac.nz) if you need any help with setting up and running this code, or with interpreting the results.
