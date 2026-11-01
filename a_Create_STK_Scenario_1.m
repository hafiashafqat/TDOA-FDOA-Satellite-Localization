%% ========================================================================
%% MODULE 1: INITIALIZE STK ENGINE & SCENARIO
%% ========================================================================
clc;
clear;
close all;

% 1. Connect to STK 11 or launch a new instance
fprintf('Connecting to STK 11 Engine...\n');
try
    app = actxGetRunningServer('STK11.Application');
catch
    app = actxserver('STK11.Application');
end

% Ensure the STK GUI is visible for visual verification
app.Visible = 1; 

% Obtain the Object Model root
root = app.Personality2;

% 2. Safe Guard: Check if a scenario is already open, and close it cleanly
if root.Children.Count > 0
    fprintf('An active scenario was detected. Unloading it to avoid conflicts...\n');
    root.CloseScenario(); 
end

% 3. Create a fresh scenario
scenarioName = 'RF_Localization';
fprintf('Creating new scenario: %s...\n', scenarioName);
root.NewScenario(scenarioName);
scenario = root.CurrentScenario;

% 4. Establish Simulation Timeline (UTC)
scenario.SetTimePeriod('06 Oct 2026 00:00:00.000', '07 Oct 2026 00:00:00.000');

% 5. Rewind animation clock to start time
root.Rewind;

% 6. Confirm success
fprintf('\n===========================================\n');
fprintf(' STK CONNECTION & INITIALIZATION SUCCESSFUL\n');
fprintf('===========================================\n');
fprintf('Active Scenario : %s\n', scenario.InstanceName);
fprintf('Simulation Start: %s\n', scenario.StartTime);
fprintf('Simulation Stop : %s\n', scenario.StopTime);
fprintf('===========================================\n');

disp('MODULE 1 COMPLETE');


%% ========================================================================
%% MODULE 1: STK SCENARIO INSPECTION
%% ========================================================================
fprintf('\nExecuting Module 2: Scenario Inspection Verification...\n');

try
    % Safely query the existing handle
    appInspection = actxGetRunningServer('STK11.Application');
    rootInspection = appInspection.Personality2;
    
    if rootInspection.Children.Count == 0
        error('Inspection failed: No active scenario loaded in the STK instance.');
    end
    
    currentScen = rootInspection.CurrentScenario;
    
    fprintf(' INSPECTION SUCCESSFUL\n');
    fprintf(' Verified Scenario Name: %s\n', currentScen.InstanceName);
    fprintf(' Verified Time Window : %s to %s\n', currentScen.StartTime, currentScen.StopTime);
    disp('MODULE 1 PASSED');

catch ME
    fprintf(' Inspection critical failure: %s\n', ME.message);
end
