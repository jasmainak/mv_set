
%Code for various classification/set estimation problems using dyadic decision
%trees (DDTs). Developed by Clayton Scott, Rice University, 2005. 

ddt_demo.m

%Demonstration script. Be sure to try different datasets on different 
%learning problems.

ddt_learn.m

%Takes training data and produces classifier or set estimate. Type "help
%ddt_learn" for options.
    
ddt_core.m

%Help file for ddt_core executable

ddt_core.dll (or other executable)

%C++ routine for minimizing an additive functional over the class of DDTs. 
%Implements algorithm described in the paper "Oracle bounds and exact 
%algorithm for dyadic classification trees," In Proceedings of the 17th. 
%Conference on Learning Theory (COLT 2004),2004. 
%Springer Lecture Notes in Artificial Intelligence (3120), 378-392, 2004. 

mvg.m

%Generates multivarate Gaussian data for testing

banana.mat

%Copy of Raetsch's synthetic banana data.