function [T,Y] = model_database(model_handle, k1, k2, k3, v1, r1, r2, a, b, t2, initial_values_nd, duration)

    [T,Y] = ode15s(@biomin_K3_IR_v3_2, [0, duration], initial_values_nd);

    
    % Modifed Hill function
    function hill = H(x)
        hill = b/(b+x^a);
    end

    % Biomineralization model structure
    function Zdot = biomin_K3_IR_v3_2(t,z) %z=(x1,x2,I,N,y)
        
        Zdot=zeros(5,1);
        Zdot(1) = -k1 * z(1);
        Zdot(2) = k1 * z(1);
        Zdot(3) = v1 * z(1) - r1 * z(2) * z(3) - t2 * (k3 * H(z(3)) * z(4)) * z(3); % dI/dt = v1.x1 - r1.x2.I - t2.dy/dt.I
        Zdot(4) = k2 * (k1*z(1)) - r2 * z(4) * (k3*H(z(3))*z(4));
        Zdot(5) = k3 * H(z(3)) * z(4);
    end


end
 