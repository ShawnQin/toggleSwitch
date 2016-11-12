# toggleSwitch
##About this program
Here is the brief description of what this program did. For the _toggle switch _ model, there is a integrated program `intrinsicExtrinsicToggleSwitch.m`. It contain all the _fucntions_ needed to do the simulation. There are three different _linear noise approximation_ functions dealing with different amplitude, correlation time scale of extrinsic fluctuation and also the system size effect. We use _Chemical Langevin Equations_ to do the simulation, although _Gillespie Algorithms_ gives qualitatively the same results

## explain of newly added functions
`toggleSwiLNAExtri.m` input of different source of extrinsic noise and return and plot the calcualted EWS based on LNA. When only one source of extrinsic noise is added, it will add two same extrinsic noise to the corresponding postion of tow SDEs. When more than one sources of extrisic noise is specified, the amlitude and correlation time scale of extrinsic noise also have to specified explicitly.
