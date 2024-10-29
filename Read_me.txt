Guide to running sims:

1. Open SCAAwAE_Sim_Functions_Sourced.R
2. Adjust your working directory path (line 7)
3. To switch between multinomial or dirichlet multinomial comment out/uncomment lines 10-11 and 92-93
4. Set number of iterations (line 13)
5. Set max number of jitters (line 14)
6. If you only want to run a limited number of scenarios you can adjust them in line 174, or comment out that line to run all scenarios.
8. Adjust number of cores to run sims on in line 181, currently set to number available - 1.
7. Run script

To check convergence:
1. Open Analysis.R
2. Adjust your working directory path (line 5)
3. To check just Hessian PD uncomment line 24, to check Hessian PD and max_grad <= 0.1 uncomment line 25
4. Run script to line 72

Plot time series of relative error:
1. Open Analysis.R
2. Adjust your working directory path (line 5)
3. Run lines 1-11
4. Run lines 79-84