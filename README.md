# GEM-A-GEneral-Memristive-transistor-model

GEM (General Memristive Transistor Model) is a general model for simulating memristive transistors. The GEM model features a time-dependent differential equation, a voltage-controlled moving window function, and a nonlinear current output function. This combination enables precise simulation of both the switching and output characteristics of memristive transistors. For more details, refer to our preprint here (https://arxiv.org/abs/2408.15140).

## Contents ##
This repository includes the GEM model for the SPICE simulator, as well as the optimization code required to tune GEM parameters based on experimental data. The folders are organized as follows:
1. Result: Contains data and figures presented in our manuscript in the Origin file format.
2. SPICE_Model: Includes the .sub file and test SPICE circuit for the GEM model.
3. Switching_Behavior: Contains data and optimization code for testing GEM’s ability to model resistance switching.
4. Physical_Limit: Includes data used to test GEM's ability to replicate changes under varying voltage stimuli amplitudes.
5. Nonlinear_Output: Contains data and optimization code for testing GEM’s nonlinear output capabilities.

## Required Software ##
1. Matlab
2. Origin
3. LTspiceXVII

## LTspice configuration instruction ##
1. Download the binary file for LTspice XVII from http://ltspice.analog.com and execute the installation program. It is recommended to install the software in the directory 'C:\Program Files\LTC\LTspiceXVII' to facilitate ease of access for subsequent simulations.
2. Add the SPICE_models to the Library Search Path in LTspice XVII manually. (Detailed instructions can be found in https://uspas.fnal.gov/materials/17NIU/LTspiceXVII%20Installation.pdf.) Notice: in the first use, please modify the .lib instructions in the .asc to guarantee successful simulation in LTspice.

## Contact ##
For inquiries, please contact:
Shengbo Wang: shengb.wang@gmail.com
Shuo Gao: shuo_gao@buaa.edu.cn