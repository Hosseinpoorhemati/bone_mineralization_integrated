function [k1_nonD, k2_nonD, k3_nonD, v1_nonD, r1_nonD, r2_nonD, a_nonD, b_nonD, t2_nonD] ...
                = non_dimensionalizer(model_handle, k1, k2, k3, v1, r1, r2, a, b, t2, characteristics)
    
    x1_c = characteristics(1);
    x2_c = characteristics(2);
    I_c = characteristics(3);
    N_c = characteristics(4);
    y_c = characteristics(5);
    t_c = 1; 

    if model_handle == "K3_IR_v1_0"
        k1_nonD = k1 * t_c ;
        v1_nonD = v1 * x1_c * t_c / I_c;
        r1_nonD = r1 * x2_c * t_c;
        r2_nonD = r2 * y_c;
        k2_nonD = k2 * x2_c / N_c; 
        k3_nonD = k3 * N_c * t_c / y_c;
        a_nonD = a;
        b_nonD = b/(I_c^a);
        t2_nonD = 1;

        % x_hat = 10^6;
        % y_hat = 10^9;
        % k1_nonD = k1;
        % v1_nonD = v1;
        % r1_nonD = r1 * x_hat;
        % r2_nonD = r2 * y_hat;
        % k2_nonD = k2; 
        % k3_nonD = k3 * (x_hat/y_hat);
        % a_nonD = a;
        % b_nonD = b/(x_hat^a);


    elseif model_handle == "K3_IR_v3_2"
        k1_nonD = k1 * t_c ;
        v1_nonD = v1 * x1_c * t_c / I_c;
        r1_nonD = r1 * x2_c * t_c;
        r2_nonD = r2 * y_c;
        k2_nonD = k2 * x2_c / N_c; 
        k3_nonD = k3 * N_c * t_c / y_c;
        a_nonD = a;
        b_nonD = b/(I_c^a);
        t2_nonD = t2 * y_c;


    end
end