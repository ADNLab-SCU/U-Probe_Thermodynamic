
% Date -> 2023-07-11

function SimY_cmTER = Fn_Fit_dGshift(dG_shift, dG_cmTER)

fbar = waitbar(0, 'Function Loaded');

R = 1.987*10^-3;  % kcal*K^-1*mol^-1
T = 273.15 + 37;  % Celcius to Kelvin
G2K = @(G) exp(-G./(R*T)); % function handle to calculate K

Exp_Conc = [250/20, 9]; % for [Target] and [Free_Q], respectively
Conc_T = Exp_Conc(1);
Conc_freeQ = Exp_Conc(2);

K_cmTER = G2K(dG_shift + dG_cmTER);

SimY_cmTER = zeros(size(dG_cmTER));

% Reaction Model
syms x
for i = 1: length(K_cmTER)
    Eqn_TER = (x).*(Conc_freeQ+x)./( (1-x).*(Conc_T-x) ) == K_cmTER(i);

    solx = double( solve(Eqn_TER, x) );
    SimY_cmTER(i) = solx(solx>0 & solx<1);

    waitbar(i/length(K_cmTER), fbar, 'Calculating Yields ~ dG_cmTER')
end
close(fbar)

end






