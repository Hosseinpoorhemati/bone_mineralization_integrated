% ISF equilibrium calculator FUNCTION

function [equilib_corr_conc, pH_init, pH_final, saturation_ratio, precipitation_rate]  = ISF_function(master_conc, iteration_value)

format compact;
format short;

% Entering the total concentrations measured in blood
pH = master_conc(1);
c(1) = master_conc(2);  %TCO3;  % HCO3 -
c(2) = master_conc(3);  %TPO4;  % HPO4 2-
c(3) = master_conc(4);  %TCa;   % Ca 2+
c(4) = master_conc(5);  %TMg;   % Mg 2+
c(5) = master_conc(6);  %TNa;   % Na +
c(6) = master_conc(7);  %TCl;   % Cl -
c(7) = master_conc(8);  %TK;    % K +
c0 = c;


f = @equations;
J = @Jacobian;

gamma_previous = ones (1,30) .* 0.5;
gamma_all = gamma_previous;

z0 = [0.001; 0.001; 0.001; 0.001; 0.001; 0.001; 0.001]; %initial guess for NR system
NR_Eps = 10^-8; % the precision of NR
act_Eps = 10^-8; % this values is for the activity calculation precision

[z, Iterations] = Newton_sys(f, J, z0, gamma_previous, c, NR_Eps, pH);

corrected = false;
counter1 = 0; % counter of how many times the following activity function runs before it finds the right values
while corrected == false
    [z, Iterations] = Newton_sys(f, J, z0, gamma_previous, c, NR_Eps, pH);
    z0 = z;
    [eq_conc, TH] = concentration (z, gamma_previous, pH);
    [gamma_all, IS] = activity(eq_conc);
    if (norm(gamma_all - gamma_previous, inf) <= act_Eps)       
        corrected = true;        
    end
    gamma_previous = gamma_all;
    counter1 = counter1 + 1;   
end

[ratio, status] = saturation(eq_conc, gamma_previous);

pH2 = -log10(eq_conc(1)*gamma_all(1)); % pH calculated from the corrected concentration
%fprintf("Initial calculated pH is %.2d \n", pH2);
pH_init = pH2;

N = iteration_value;
plot_ratio = zeros(1,N);
plot_IS = zeros(1,N);
plot_gamma = zeros(N,length(eq_conc));
plot_status = zeros(1,N);
plot_TH = zeros(1,N);
plot_C = zeros(N,7);
plot_conc = zeros(N,length(eq_conc));
plot_pH = zeros(1,N);
rates = zeros (1,N);
times = zeros (1,N);

S_threshhold = 1; 

counter2 = 1; % the counter used to loop the following while section
while counter2<= N && ratio > S_threshhold

    HAP_ion = [eq_conc(10), eq_conc(9), eq_conc(2)]; % vector of current Ca, PO4, OH conc.  
    [live_rate, change] = new_precipitation_1(HAP_ion, gamma_previous);
    rates(counter2) = live_rate(2);
    times(counter2) = live_rate(1);

    % next we must update the 4 mass balance  
    if ratio > S_threshhold            
        c(3) = c(3) - change(1); % Ca2 +
        c(2) = c(2) - change(2); % HPO4 2-
    elseif ratio <= S_threshhold
        fprintf("Less or equal to threshold - Break!")
        break
        %c(3) = c(3) + change(1); % Ca2 +
        %c(2) = c(2) + change(2); % HPO4 2-  
    end

   if c(2)< 0 
        fprintf('\n P depleted! \n')
        counter2;
        break
    end

    corrected = false;
    counter3 = 0; % counter of how many times the following activity function runs before it finds the right values
    while corrected == false
        [z, Iterations] = Newton_sys(f, J, z0, gamma_previous, c, NR_Eps, pH);
        z0 = z;
        [eq_conc, TH] = concentration (z, gamma_previous, pH);
        [gamma_all, IS] = activity(eq_conc);
        if (norm(gamma_all - gamma_previous, inf) <= act_Eps) %|| (activity_correction_iter > 4)        
            corrected = true;        
        end
        gamma_previous = gamma_all;
        counter3 = counter3 + 1;
    end
    
    [ratio, status] = saturation(eq_conc, gamma_previous);
    pH2 = -log10(eq_conc(1)*gamma_all(1));
  
    plot_pH(counter2) = pH2;
    plot_ratio(counter2) = ratio;
    plot_IS(counter2) = IS;
    plot_gamma(counter2,:) = gamma_previous;
    plot_status(counter2) = status;
    plot_TH(counter2) = TH;
    plot_C(counter2, :) = c;
    plot_C_normalized(counter2, :) = c./c0; % total concentration of major components available in solution divided by their initial blood level
    plot_conc(counter2, :) = eq_conc;
    counter2 = counter2 + 1;
end

equilib_corr_conc = plot_conc;
pH_final = pH2;
saturation_ratio = ratio;
precipitation_rate = rates(1);

% if iteration_value > 1
% 
%     %%%%%%% FIGURE 6-A of MANUSCRIPT %%%%%
%     figure(11)
%     hold on
%     times_ax = zeros(1,length(times));
%     times_ax(1) = times(1);
%     for j = 2:length(times)
%         times_ax(j) = times(j) + times_ax(j-1);
%     end
%     times_ax = times_ax / (60*60); % to be in hours
% 
%     plot(times_ax(1:length(times)), rates(1:length(times)))
%     legend('pH 7.3','pH 7.4','pH 7.55')
%     ylabel("HAP Precipitation Rate (Mol/(L.s))")
%     xlabel("Time (h)")
% 
%     %%%%%%% FIGURE 6-B of MANUSCRIPT %%%%%
%     figure(12)
%     hold on
%     % in a 1 by 1 by 1 um3 cube:
%     volume = (10^-6)^3 * 1000; % unit: litre
%     molar_weight = 502.3; % g/mol
%     rate_in_gram_per_second = rates .* (volume * molar_weight);
%     accumulative_mass = zeros (1,length(rate_in_gram_per_second));
%     accumulative_mass(1) = rate_in_gram_per_second(1);
%     for m = 2:length(rate_in_gram_per_second)
%         accumulative_mass(m) = rate_in_gram_per_second(m) + accumulative_mass(m-1);
%     end
% 
%     plot(times_ax, accumulative_mass);
%     legend('pH 7.3','pH 7.4','pH 7.55')
%     xlabel('Time (h)')
%     ylabel("HAP Mass (gr)")
% 
% 
%     figure(17)
%     plot (1:counter2-1,plot_conc(1:counter2-1,10), '--')
%     ylabel("iCa")    
% 
% end 


function [z, Iters] = Newton_sys(f, J, z0, g, c, Eps, pH)
    iJ = @(z) inv(J(z,g,c,pH));
    z = z0 - iJ(z0)*f(z0,g,c,pH) ;
    
    for i=1:length(c)
        if z(i)<0
            z(i) = z0(i)/10;
        end 
    end
    
    Iters = 0 ;
    while norm(z-z0, inf) >= Eps
        z0 = z ;
        z = z0 - iJ(z0)*f(z0,g,c,pH) ;
        
        for i=1:length(c)
            if z(i)<0
                z(i) = z0(i)/10;
            end 
        end
        
        Iters = Iters + 1 ;
    end
end


function J = Jacobian(x,g,c,pH)
    J = zeros (length(x), length(x));
    for i = 1:length(x)
        dx = x;
        eps = 10^-8 * x(i);
        dx(i) = dx(i) + eps; 
        J(:,i) = (equations(dx,g,c,pH) - equations(x,g,c,pH)) / eps;
    end
end 
 

function F = equations (z, g, c, pH)
    % g is the activity coeeficient vector
    % C is the blood concentration of measured components
    
    % Defining equilibrium values
    K(1) = 10^-14;                      % water
    K(2) = 10^-6.31;                    % H2CO3
    K(3) = 10^-10.25;                   % HCO3 - 
    K(4) = 10^-2.196;                   % H3PO4
    K(5) = 10^-7.185;                   % H2PO4 2
    K(6) = 10^-12.19;                   % HPO4 3-
    K(7) = 10^1.16;                     % CaHCO3 +
    K(8) = 10^3.38;                     % CaCO3
    K(9) = 25.12;                       % CaOH +
    K(10) = 31.9;                       % CaH2PO4 +
    K(11) = 6.81*10^2;                  % CaHPO4
    K(12) = 3.46*10^6;                  % CaPO4 -
    K(13) = 10^0.62;                    % MgHCO3 +
    K(14) = 10^1.87;                    % MgCO3
    K(15) = 10^2.19;                    % MgOH +
    K(16) = 10^0.4;                     % MgH2PO4 +
    K(17) = 10^1.8;                     % MgHPO4
    K(18) = 10^3.3;                     % MgPO4 -
    K(19) = 0.21;                       % NAHPO4 -
    K(20) = 10^-6.82;                   % NaH2PO4
    K(21) = 1/29.3;                     % NaCl
    K(22) = 2.5;                        % KHPO4 -
    H_ion = 10^-pH;
    
    % y(5) => z(1)
    % y(8) => z(2)
    % y(10) => z(3)
    % y(17) => z(4)
    % y(24) => z(5)
    % y(27) => z(6)
    % y(29) => z(7)
    
    H_ion = H_ion*g(1); 
    z(1) = z(1)*g(5);
    z(2) = z(2)*g(8);
    z(3) = z(3)*g(10);
    z(4) = z(4)*g(17);
    z(5) = z(5)*g(24);
    z(6) = z(6)*g(27);
    z(7) = z(7)*g(29);
      
    F(1) = (K(3)/K(2)) * H_ion * z(1) / g(3) + (K(3)*z(1)/H_ion)  / g(4) + z(1) / g(5) ...
        + (K(7)*z(1)*z(3)) / g(11)  + (K(3)*K(8)*z(1)*z(3)/H_ion) / g(12)  + ...
        (K(13)*z(1)*z(4)) / g(18)  + (K(3)*K(14)*z(1)*z(4)/H_ion) / g(19)  - c(1);
    
    F(2) = ((1/(K(4)*K(5)))*((H_ion)^2)*z(2)) / g(6)  + ((1/K(5))*H_ion*z(2)) / g(7)...
        + z(2) / g(8)  + (K(6)*z(2)/H_ion) / g(9)  + ((K(10)/K(5))*H_ion*z(2)*z(3)) / g(14) ...
        + (K(11)*z(2)*z(3)) / g(15) + (K(6)*K(12)*z(2)*z(3)/H_ion) / g(16)  + ...
        ((K(16)/K(5))*H_ion*z(2)*z(4)) / g(21)  + (K(17)*z(2)*z(4)) / g(22)  + ...
        (K(6)*K(18)*z(2)*z(4)/H_ion) / g(23)  + (K(19)*z(2)*z(5)) / g(25)  ...
        + ((K(20)/K(5))*H_ion*z(1)*z(5)) / g(26)  + (K(22)*z(2)*z(7)) / g(30)  - c(2);
    
    F(3) = z(3) / g(10)  + (K(7)*z(1)*z(3)) / g(11)  + (K(3)*K(8)*z(1)*z(3)/H_ion) / g(12)  + ...
        (K(1)*K(9)*z(3)/H_ion) / g(13)  + ((K(10)/K(5))*H_ion*z(2)*z(3)) / g(14)  + ...
        (K(11)*z(2)*z(3)) / g(15)  + (K(6)*K(12)*z(2)*z(3)/H_ion) / g(16)  - c(3);
    
    F(4) = z(4) / g(17)  + (K(13)*z(1)*z(4)) / g(18)  + (K(3)*K(14)*z(1)*z(4)/H_ion) / g(19)  + ...
        (K(1)*K(15)*z(4)/H_ion) / g(20)  + ((K(16)/K(5))*H_ion*z(2)*z(4)) / g(21)  + ...
        (K(17)*z(2)*z(4)) / g(22)  + (K(6)*K(18)*z(2)*z(4)/H_ion) / g(23)  - c(4);
    
    F(5) = z(5) / g(24)  + (K(19)*z(2)*z(5)) / g(25)  + ((K(20)/K(5))*H_ion*z(1)*z(5)) / g(26) ...
        + (K(21)*z(5)*z(6)) / g(28)  - c(5);
    
    F(6) = z(6) / g(27)  + (K(21)*z(5)*z(6)) / g(28)  - c(6);
    
    F(7) = z(7) / g(29)  + (K(22)*z(2)*z(7)) / g(30)  - c(7);    
    
    F = F.';
end


function [eq_conc, TH] = concentration (z, gamma, pH)
    % caclulte concentration of each component using the output of Newton_sys
    % function and the activities. input arguments are:
    % x = concentrations of 7 main components
    % gamma = gamma for all components
    % output: uncorrected but equilibrium concentrations for all
    
    % Defining equilibrium values
    K(1) = 10^-14;                      % water
    K(2) = 10^-6.31;                    % H2CO3
    K(3) = 10^-10.25;                   % HCO3 - 
    K(4) = 10^-2.196;                   % H3PO4
    K(5) = 10^-7.185;                   % H2PO4 2
    K(6) = 10^-12.19;                   % HPO4 3-
    K(7) = 10^1.16;                     % CaHCO3 +
    K(8) = 10^3.38;                     % CaCO3
    K(9) = 25.12;                       % CaOH +
    K(10) = 31.9;                       % CaH2PO4 +
    K(11) = 6.81*10^2;                  % CaHPO4
    K(12) = 3.46*10^6;                  % CaPO4 -
    K(13) = 10^0.62;                    % MgHCO3 +
    K(14) = 10^1.87;                    % MgCO3
    K(15) = 10^2.19;                    % MgOH +
    K(16) = 10^0.4;                     % MgH2PO4 +
    K(17) = 10^1.8;                     % MgHPO4
    K(18) = 10^3.3;                     % MgPO4 -
    K(19) = 0.21;                       % NAHPO4 -
    K(20) = 10^-6.82;                   % NaH2PO4
    K(21) = 1/29.3;                     % NaCl
    K(22) = 2.5;                        % KHPO4 -
    H_ion = 10^-pH;
      
    y(1) = H_ion*gamma(1); 
    y(5) = z(1)*gamma(5);
    y(8) = z(2)*gamma(8);
    y(10) = z(3)*gamma(10);
    y(17) = z(4)*gamma(17);
    y(24) = z(5)*gamma(24);
    y(27) = z(6)*gamma(27);
    y(29) = z(7)*gamma(29);
    
    
    y(2) = (K(1)/y(1)) / gamma(2);
    y(3) = ((K(3)/K(2)) * y(1) * y(5)) / gamma(3) ;
    y(4) = (K(3)*y(5)/y(1)) / gamma(4);
    y(6) = ((1/(K(4)*K(5)))*(y(1)^2)*y(8)) / gamma(6);
    y(7) = ((1/K(5))*y(1)*y(8)) / gamma(7);
    y(9) = (K(6)*y(8)/y(1)) / gamma(9);
    y(11) = (K(7)*y(5)*y(10)) / gamma(11);
    y(12) = (K(3)*K(8)*y(5)*y(10)/y(1)) / gamma(12);
    y(13) = (K(1)*K(9)*y(10)/y(1)) / gamma(13);
    y(14) = ((K(10)/K(5))*y(1)*y(8)*y(10)) / gamma(14);
    y(15) = (K(11)*y(8)*y(10)) / gamma(15);
    y(16) = (K(6)*K(12)*y(8)*y(10)/y(1)) / gamma(16);
    y(18) = (K(13)*y(5)*y(17)) / gamma(18);
    y(19) = (K(3)*K(14)*y(5)*y(17)/y(1)) / gamma(19);
    y(20) = (K(1)*K(15)*y(17)/y(1)) / gamma(20);
    y(21) = ((K(16)/K(5))*y(1)*y(8)*y(17)) / gamma(21);
    y(22) = (K(17)*y(8)*y(17)) / gamma(22);
    y(23) = (K(6)*K(18)*y(8)*y(17)/y(1)) / gamma(23);
    y(25) = (K(19)*y(8)*y(24)) / gamma(25);
    y(26) = ((K(20)/K(5))*y(1)*y(5)*y(24)) / gamma(26);
    y(28) = (K(21)*y(24)*y(27)) / gamma(28);
    y(30) = (K(22)*y(8)*y(29)) / gamma(30);
    
    y(1) = y(1)/gamma(1); 
    y(5) = y(5)/gamma(5);
    y(8) = y(8)/gamma(8);
    y(10) = y(10)/gamma(10);
    y(17) = y(17)/gamma(17);
    y(24) = y(24)/gamma(24);
    y(27) = y(27)/gamma(27);
    y(29) = y(29)/gamma(29);

    eq_conc = y;
    
    TH = y(1) + y(2) + 2*y(3) + y(5) + 3*y(6) + 2*y(7) + y(8) + y(11) + y(13) + ...
        2*y(14) + y(15) + y(18) + y(20) + 2*y(21) + y(22) + y(25) + 2*y(26) + y(30);
end


function [gamma, IS] = activity(concentrations)            
    z = [1, -1, 0, -2, -1, 0, -1, -2, -3, 2, 1, 0, 1, 1, 0, -1, 2, 1, 0, 1, 1, 0, -1, 1, -1, 0, -1, 0, 1,-1]; 
    IS = 0;        
    for i = 1:30
        IS = IS + 0.5 * (concentrations(i)*(z(i)^2));
    end

    temp = 37;
    A = 0.486 + (6.07 * (10^-4) * temp) + (6.43 * (10^-6) * (temp^2));      
    gamma = zeros(1,30);
    for i= 1:30
        gamma(i) = 10 ^ (-A * (z(i)^2) * ((sqrt(IS)/(1+sqrt(IS)))-(0.3*IS)));
    end       
end


function [sat_ratio, sat_status] = saturation(y, g)
    % calculates the saturation ratio and provides a True/False output
    
    % The precipitation component(s): HAP (hydroxyapatite): Ca10 (PO4)6 (OH)2 or Ca5 (PO4)3 OH
    % Reaction: 5Ca2+ + 3PO43- + OH-  =>  Ca5(PO4)3OH      with KSP = 2.35 * 10^-59
    
    % The concentrations are corrected and come from the concentration function output and:
    % Ca = y(10)
    % OH = y(2)
    % PO4 = y(9)
    
    % S method 1 - considers the saturation based on HAP
    IP = ((g(10)*y(10))^5) * ((g(9)*y(9))^3) * g(2)*y(2);
    Ksp = 2.03 * 10^-59; % from paper I-39-Sohnet 2011 @ 37 C
    sat_ratio = (IP/Ksp)^(1/9); % Intoduced in paper I-39: Sohnet 2011: this is in fact solution supersaturation, expected to be in the order of 12-13 
    %sat_ratio = IP/Ksp; % this is what we cal saturation ratio, >1 means supersaturated. 
    %sat_ratio = log10(IP/Ksp); % this is the saturation index, >0 means supersaturated

    %S method 2 - considers the saturation based on ACP not HAP
    %IP = ((g(10)*y(10))^3) * ((g(9)*y(9))^2);
    %Ksp = 10^-26; 
    %sat_ratio = (IP/Ksp)^(1/5);

    if sat_ratio > 1
        sat_status = +1;
    elseif sat_ratio < 1
        sat_status = -1;
    else
        sat_status = 0;
    end
end


function [live_rate, conc_change] = new_precipitation_1(HAP_ion, g)
    % This mode uses the R = Ks*s*gamma2*gamma3*[Ca]*[PO4] to calculate
    % precipitation
    Ca_act = HAP_ion(1) * g(10);
    P_act = HAP_ion(2) * g(9);

    k = 173.4; % rate constant
    s = 14.21; % m^2 / g
    rate = k * s * Ca_act * P_act; % mol HAP / L.s = molar/s = concentration per s
   
    k1 = 10^-6;
    s1 = 50; 
    OH_act = HAP_ion(3) * g(2);
    Ksp = 2.03 * 10^-59;
    rate_k = k1*s1* (((Ca_act^5)*(P_act^3)*(OH_act))^(1/9) - Ksp^(1/9))^1.25;
    
    t = 1  ; % time step, equal to 1 sec
    delta_HAP = rate * t; % molar change in HAP in time t 

    live_rate = [t,rate]; % unit = second
    conc_change = [5*delta_HAP, 3*delta_HAP, delta_HAP]
end


end 
