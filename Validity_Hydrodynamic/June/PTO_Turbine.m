function [P_p, P_t, T_t, eta,m_t,valve,Psi,Phi,Tau] = PTO_Turbine(p_star,Omega,u_v,Int,WEC,Turb,CL_num)

Omega = real(Omega);

% Pressure across the turbine
p_turbine = p_star * WEC.p_at;

% Turbine characteristic curves
Psi = abs(  p_turbine / (WEC.rho_at * (Turb.D)^2 * Omega^2 ));

switch Turb.Turb 
    % Wells
    case 1
        valve = PTO_valve(Omega,u_v,CL_num,Int,Turb);
        radius = 1;
        
%         Psi = max(1e-4, abs(p_turbine) / (WEC.rho_at * (Turb.D*radius)^2 * Omega^2 ));
%         if Psi > 2.5, Psi = 2.5; end
        
        %Flow rates       
        Phi = valve * Turb.Kt * Psi; 
        if isnan(Phi), Phi = 0; end
        
        % Flow rates        
        Q_t = valve * sign( p_star ) * Omega *  (Turb.D*radius)^3 * Phi;        
        m_t = WEC.rho_at * Q_t;    
        
         Tau = interp1(Turb.Phi_tbl,Turb.Tau_tbl,Phi,'linear');
%         if isnan(Tau), Tau = 0; end
%         
%         T_t = Tau * WEC.rho_at * (Turb.D*radius)^5 * Omega^2;
%         P_t = T_t * Omega;
% 
%         eta = P_t/(p_turbine * Q_t);
%         eta = Tau  / (Psi * Phi);
        eta  = interp1(Turb.Phi_tbl,Turb.Eta_tbl,Phi,'linear');
        if eta < -1, eta=-1; elseif isnan(eta), eta = 0; end       
        P_t = eta * p_turbine * Q_t;
        T_t = P_t / Omega;
        
    % Axial impulse 
    case 2  
        P_t = 0; T_t = 0; eta = 0; m_t = 0; valve = 1; Tau = 0;

    % Biradial
    case 3
        

        if( Psi > Turb.Psi_max)
            Psi = Turb.Psi_max; 
        end
        
        % Flow rates
        valve = PTO_valve(Omega,u_v,CL_num,Int,Turb);        
        Phi = valve * sign(Psi)*0.195 * (abs(Psi)+0.000001)^0.6;
        
        Q_t = valve * sign( p_star ) * Omega *  Turb.D^3 * Phi;        
        m_t = WEC.rho_at * Q_t;       
        
        if( Psi > 3.3116337496519517 )
            eta = 1.0/(1.442750893460873 + 0.2066910461657107*Psi); 
        elseif( Psi > 0.35075016168255563 )
            idx = 1;
            eta = Turb.SP(idx,4) + Turb.SP(idx,3)*Psi + Turb.SP(idx,2)*Psi^2 + Turb.SP(idx,1)*Psi^3;
        elseif( Psi > 0.3082881279543439 )
            idx = 2;
            eta = Turb.SP(idx,4) + Turb.SP(idx,3)*Psi + Turb.SP(idx,2)*Psi^2 + Turb.SP(idx,1)*Psi^3;
        elseif( Psi > 0.18131721834023762 )
            idx = 3;
            eta = Turb.SP(idx,4) + Turb.SP(idx,3)*Psi + Turb.SP(idx,2)*Psi^2 + Turb.SP(idx,1)*Psi^3;
        else
            idx = 4;
            eta = Turb.SP(idx,4) + Turb.SP(idx,3)*Psi + Turb.SP(idx,2)*Psi^2 + Turb.SP(idx,1)*Psi^3;
        end               

        % Turbine power and torque
        P_t = eta *  WEC.p_at * p_star * Q_t;
        T_t = P_t / Omega; 
        Tau=0;
        
    otherwise 
        P_t = 0; T_t = 0; eta = 0; m_t = 0; valve = 1;Tau=0;
end
P_p = p_turbine * Q_t;
end

function v_pos = PTO_valve(Omega,u_v,CL_num,Integ,Turb)

if Integ == 1, v_pos = u_v;
else
    % Safety valve closes when maximal turbine speed reached ...
    % ... until speed less than a cut in speed
    if u_v == 0 && Omega <= Turb.Om_ci
        v_pos = 1;  
    elseif Omega > Turb.Om_co || u_v == 0   
        v_pos = 0;     
    else v_pos = u_v;
    end
end

if CL_num == 0, v_pos = 1; end

end