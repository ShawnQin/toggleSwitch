# toggleSwitch
##About this program
Here is the brief description of what this program did. For the _toggle switch _ model, there is a integrated program `intrinsicExtrinsicToggleSwitch.m`. It contain all the _fucntions_ needed to do the simulation. There are three different _linear noise approximation_ functions dealing with different amplitude, correlation time scale of extrinsic fluctuation and also the system size effect. We use _Chemical Langevin Equations_ to do the simulation, although _Gillespie Algorithms_ gives qualitatively the same results

## explain of newly added functions in toggle switch model
`toggleSwiLNAExtri.m` input of different source of extrinsic noise and return and plot the calcualted EWS based on LNA. When only one source of extrinsic noise is added, it will add two same extrinsic noise to the corresponding postion of tow SDEs. When more than one sources of extrisic noise is specified, the amlitude and correlation time scale of extrinsic noise also have to specified explicitly.

`ExportDataR.m`  integrates all the stochastic simulation results from the CLS servers, out put a neat data struct that can be further used in R to plot figures. Two data structs are generated _ExtriFig1SimuAmp.mat_ and _ExtriFig2SimuTime.mat_ which represents different amplitudes and correlation time scales.

`toggSwiCLEStaticEnsembleExt.m` this program tries to demonstrate the EWS got from static ensemble data. Basically, after simulating an ensemble of repeats and only take serveral data point of the trajectory to do all the statistics.

`toggTimeDetrendExtri.m` this program the detrending of single time serials won't eliminate the effects of extrinsic noise. When control parameter changes gradully from one regime to another, the trend of EWS is quantified by Kandell's tau correaltion. Since I also used this program to test Chen's model of large noise, so it might look confusing.

## functions in negative feedback model
`negativeVarCorrCLEext.m` this program implement CLE algorithm of a negative feedback motif and shows the false positve resutls when extrinsic noise is considered.



