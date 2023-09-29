%% This script calls getPlots.m to plot and save the phase portrait, saddle connector solution, probability density, and parametric plots

clear all
close all
clc
%% Inputs for the initial s-value, winding number, energy guess, gamma value, and folder name

s0 = 0; % initial value of s
n = 3;
E = 0.9; % value of constant E
a = 7.5; % value of constant a
FolderName = [.];   % Your destination folder


figures = getPlots(n,s0,E,a);
%% Saves figures as .png to a specified folder
FigList = findobj('Type', 'figure');
Name_list = ["parametric" "prob" "saddle" "phaseportrait"];

for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = ["gamma_" + a + "_" + Name_list(iFig)];
  saveas(FigHandle, FolderName + "\" + FigName + "_WN" + n + ".png");
end