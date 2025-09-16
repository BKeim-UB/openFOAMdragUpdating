# openFOAMdragUpdating


This repository is used to update two aspects of particle-particle drag models in OpenFOAM.

The DragModels folder contains a carnahan-starling formulation for particle-particle drag, which can be added to  /applications/OpenPDAC/phaseSystem/interfacialModels/dragModels. These files contain a 'baked in' radial distribution function that is calculated based on the two phases.

Each indivudial phase, and the granular pressure and viscosity, also use a radial distribution function which is located at applications/OpenPDAC/momentumTransportModels/kineticTheoryModels/radialModel. The folder here contains another carnahan-starling radial distribution function so that each individual phase, and the interaction between phases uses the same style of radial distribution function, which is mathematically self-similar (e.g., the polydisperse case collapses to the monodisperse case when the number of phases N=1.

