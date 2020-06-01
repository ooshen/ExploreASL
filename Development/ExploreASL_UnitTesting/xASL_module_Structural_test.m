% xASL_module_Structural_test Script to test the xASL_module_Structural function
%
% FORMAT:           RESULT = runtests('xASL_module_Structural_test');
%
% INPUT:            None
%
% OUTPUT:           Console window log
%
% OUTPUT FILES:     None
%
% -----------------------------------------------------------------------------------------------------------------------------------------------------
% DESCRIPTION: This script is supposed to run various unit tests of xASL_module_Structural:
%
%           1) Run a test using the default TestDataSet inputs with high quality setting
%           2) Run a test with ...
%
% EXAMPLE: RESULT = runtests('xASL_module_Structural_test');
% __________________________________
% Copyright (C) 2015-2020 ExploreASL

% DEFAULTS

% Get ExploreASL directory
xASLdir = uigetdir(pwd,'Select ExploreASL directory...');

% Define test directory
testDir = uigetdir(pwd,'Select test directory...');
fprintf('Creating test folder in %s...\n', testDir)
mkdir(fullfile(testDir,'TestFolder'))

% Copy TestDataSet to this directory
xASL_Copy(fullfile(xASLdir,'External\TestDataSet'), fullfile(testDir,'TestFolder','TestDataSet'))
fprintf('Copy test data to %s...\n', fullfile(testDir,'TestFolder','TestDataSet'))

% PRECONDITIONS

%% Test 1: Default TestDataSet with high quality setting

% Load test inputs
testInput = load('Development\ExploreASL_UnitTesting\TestData\xASL_module_Structural_testData_Input.mat');
fprintf('Loading test input...\n')

% Modify path variables
fprintf('Modifying test input...\n')
testInput.x.MyPath = xASLdir;
testInput.x.SpaghettiDir = fullfile(testDir,'TestFolder','TestDataSet','Population','SpaghettiPlots');
testInput.x.HistogramDir = fullfile(testDir,'TestFolder','TestDataSet','Population','Histograms');
testInput.x.StatsMaps = fullfile(testDir,'TestFolder','TestDataSet','Population','StatsMaps');
testInput.x.SPMDIR = fullfile(xASLdir,'External','SPMmodified');
testInput.x.SPMpath = fullfile(xASLdir,'External','SPMmodified');
testInput.x.LockDir = fullfile(testDir,'TestFolder','TestDataSet','lock','xASL_module_Structural','Sub-001');
testInput.x.SUBJECTDIR = fullfile(testDir,'TestFolder','TestDataSet','Sub-001');

% Run test
% [result, x] = xASL_module_Structural(testInput.x);

% Use assert for outputs
% assert(isfield(x,'ModuleName'))                 % Check if new ModuleName field was created
% assert(isfield(x,'result'))                     % ...

% WORK IN PROGRESS
% Which fields always have to be in the resulting x structure?
% 
% Run again and also save 'result' variable to assert this too
% assert(isfield(result,''))  %










