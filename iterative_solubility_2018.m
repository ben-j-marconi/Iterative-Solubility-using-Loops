% Solving for gypsum solubility iteratively
clear all
close all

% In this lab, we will solve for gypsum solubility without assuming
% ideality. We will learn to use 'for' and 'while' loops in the process and
% think about how to minimize an error parameter to maximize program
% efficiency. DO NOT CHANGE THE NAMES OF SPECIFIED VARIABLES.

% The reaction under consideration is dissolution of barite in water

%% 1. Input Ksp to solve for Ksp of rxn
Ksp = 4.32*10^-5;

%% 2. Solve for the molality (Ca2+ used in ex)
% Ideal conditions assumption
% Activity coefficients for Ca and SO4
gammaCa = 1;
gammaSO4 = gammaCa;

% Concentrations (molalities) for Ca and SO4 (note these are also the
% activities for our assumptions)
mCa = (Ksp)^(1/2);
mSO4 = (Ksp)^(1/2);

%% 3. Now, let's account for non-ideal conditions.
% Calculate the activity coefficient for Ca2+ and SO42-
% Calculate ionic strength of solution
I = .5*((mCa*2^2) + (mSO4*2^2));

% Calculate activity coefficients for each ion. 
A = 0.5085;
B = 0.3281*10^8;
ao = 5*10^-8;
gammaCa = 10^((-A*2^2*sqrt(I))/(1+B*ao*sqrt(I)));
gammaSO4 = 10^((-A*2^2*sqrt(I))/(1+B*ao*sqrt(I)));

%% 4. We can now compute more accurate concentrations and ionic strength. Do so.
% Note that these are 'more accurate' - the final values must be calculated
% iteratively. The 'new_mCa' names will come into utility in a moment
mCa = (Ksp/(gammaCa*gammaSO4))^(1/2);
mSO4 = (Ksp/(gammaCa*gammaSO4))^(1/2);

%% Now, let's try the same thing, but use loops
% The general premise here is that at some point, the changes in ionic
% strength become meaningless - we have essententially found equilibrium
% conditions for all parameters. The real question is how do we know when
% to stop looping? We'll explore this using 'for' and 'while' loops.

%% 5. Use a 'for' loop to iteratively solve for solubility
% Start from ideal conditions
gammaCa = 1;
gammaSO4 = 1;
mCa = (Ksp/(gammaCa*gammaSO4))^(1/2);
mSO4 = (Ksp/(gammaCa*gammaSO4))^(1/2);

% Set up a matrix to keep track of relevant parameters
% The -9999 values are stand-ins. They are just flags for us to note that
% for this first iteration, we have not accounted for the ionic strength or
% the change in ionic strength between iterations, which is our error
% metric (note that this cannot exist for the first iteration).
% itProg_forLoop(1,:) = [I change1 mCa gammaCa mSO4 gammaSO4];
itProg_forLoop(1,:) = [-9999 -9999 mCa gammaCa mSO4 gammaSO4];

% We have technically already done iteration 1, so we'll start on 2.
counter1 = 1;
for k = 2:11 % This determines how many iterations we do (it is a known number)
    % Update counter
    counter1 = counter1 + 1;
    
    % Calculate new ionic strength of solution
    I = .5*((mCa*2^2) + (mSO4*2^2));
    
    % Now perform the iteration - note we are overwriting 'old' iterations
    % with 'new' iteration information each time!
    % Calculate new activity coefficients
    A = 0.5085;
    B = 0.3281*10^8;
    ao = 5*10^-8;
    gammaCa = 10^((-A*2^2*sqrt(I))/(1+B*ao*sqrt(I)));
    gammaSO4 = 10^((-A*2^2*sqrt(I))/(1+B*ao*sqrt(I)));
    
    % Now compute more accurate concentrations
    mCa = (Ksp/(gammaCa*gammaSO4))^(1/2);
    mSO4 = (Ksp/(gammaCa*gammaSO4))^(1/2);

    
    % Calculate the error parameter. Note that this is (new I)-(old I)
    change1 = abs(I-itProg_forLoop(counter1-1,1));
    
    % Update the dummy variable to watch iteration progress
    itProg_forLoop(counter1,:) = [I change1 mCa gammaCa mSO4 gammaSO4];
end

%% 5. Use a 'while' loop to iteratively solve for solubility
% Start from ideal conditions
gammaCa = ;
gammaSO4 = ;
mCa = ;
mSO4 = ;

% Keep track of some info for later
% The -9999 values are stand-ins. They are just flags for us to note that
% for this first iteration, we have not accounted for the ionic strength or
% the change in ionic strength between iterations, which is our error
% metric (note that this cannot exist for the first iteration).
% itProg_whileLoop(1,:) = [I change2 mCa gammaCa mSO4 gammaSO4];
itProg_whileLoop(1,:) = [-9999 -9999 mCa gammaCa mSO4 gammaSO4];

% We have technically already done iteration 1, so we'll start on 2.
counter2 = 1; % counts the number of loops

% When should the while loop stop? Set the value here - keep in mind that
% this records the change between ionic strengths of subsequent iterations
% (another reason to set the first 'I' to -9999.
minchange = ;

% The while loop will run until we reach our error metric. To start the
% loop, we have to be outside of this error (i.e., the 9999 is a
% placeholder).
change2 = 9999;

while change2 > minchange % This determines how many iterations we do (it is an unknown number)
    % Update counter
    counter2 = counter2+1;
    
    % What needs to go it in here?
    
    % Calculate the error parameter. We must be below minchange before the
    % loop will stop. Note that this is (new I)-(old I)
    change2 = abs(I-itProg_whileLoop(counter2-1,1));
    
    % Update the dummy variable to watch iteration progress
    itProg_whileLoop(counter2,:) = [I change2 mCa gammaCa mSO4 gammaSO4];
end

