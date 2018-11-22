%% Path Management
% Usage: the current folder must be the root folder. All subfolders are
% included in the path.
addpath(genpath(pwd));
%% Plot properties
% Change default plot properties
set(0,'defaultlinelinewidth',1.3); % LineWidth in point
%set(0, 'defaultTextInterpreter', 'latex');
% Dock figure by default
set(0,'DefaultFigureWindowStyle','docked');
%set(0,'DefaultFigureWindowStyle','default');
