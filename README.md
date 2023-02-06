# Heat Transfer Projects

## English

Repository with coding projects for the Heat Transfer course at Universidad de los Andes. This repository has been modified by the author to reflect some further extensions and modifications that were posteriorly performed.

The repository encompasses three projects:

1. **Project 1 - Automated design of a fin:** This project consists of a simulation of the heat profile and efficiency of a heat dissipating fin, which has some given geometric constraints. The program uses a gradient descent method as a means to maximize fin efficiency and then plots the evolution of this optimization process.

This project should receive user defined values for the physical parameters, but as of now these are hard-coded (but can be easily changed in code without affecting the model). The [optimization evolution](optimization_evolution.png) file only presents a few initial iterations of optimization, as the program has not been run completely since the recent correction of a series of bugs coming from the original version.

The program utilizes a computationally intensive method for surface area calculation which affects performance significantly; nonetheless, this method has been left in the code as it has been the only attempted strategy which provides consistently accurate results. This aspect will be re-evaluated in future versions. A bug that allows for impossible temperature profiles is present, but is only effectful with some fin profiles that possess sharp changes in cross-section. This will be fixed in a future version. An adiabatic boundary condition at the fin tip is used; changes will be made to allow for different boundary conditions to be applied. Some further changes will be made regarding the program's graphical outputs.

*A document with the mathematical formulation of the discretizations used (these are not trivial) will be uploaded soon.*

2. **Project 2 - Simulation of the heat profile in a mechanical component:** *Awaiting updated version*
3. **Project 3 - Simulation aided design of a refrigeration system:** *Awaiting updated version*

## Español

Repositorio con los proyectos para el curso de Transferencia de Calor en la Universidad de los Andes.
