% SynDiag Validation - Part 1
%
% This script should be ran everytime a function in SynDiag is modified. 
% It aims at evaluating all the diagnostics in three situations:
% -Ideal
% -Noise proportional
% -Noise proportional + Noise absolute
%
% A first check can be done just by comparing the output (figure 1) with
% the figure 1 in the document Validation_checks.docx (or
% Validation_checks.pdf)
%
% Final validations are done by SimPla module responsibles: 
% Riccardo Rossi        (r.rossi@ing.uniroma2.it)
% Simone Kaldas         (simone.kaldas@students.uniroma2.eu)
% Ivan Wyss             (ivan.wyss@uniroma2.it)
% Novella Rutigliano    (novella.rutigliano@alumni.uniroma2.eu)

%%
clear; clc;

