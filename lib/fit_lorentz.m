%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fit the data 'data' to a lorentzian:
%
%       L(x) = C - A * B^2/( (x-x0)^2 + B^2 )
%
%       BUT: The maximum (min_x,min_y) is fixed !
%
%       min_x = x-value of the lorentzian peak
%       min_y = y-value of the lorentzian peak
% optional:
%       lambda_min = lower boundary for the fit [nm]
%       lambda_max = upper boundary for the fit [nm]
%       plot_lower = lower boundary for the plot [nm]
%       plot_upper = upper boundary for the plot [nm]
%
%       FWHM = 2*B
%       Q = x0/FWHM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [A,B,C,x0,SSE] = fit_lorentz(data,min_x,min_y,B0,lambda_min,lambda_max,plot_lower,plot_upper)
    
    if nargin==4
        [A,B,C,x0,SSE] = fit_lorentz1(data,min_x,min_y,B0);
    else
        if nargin==8
            [A,B,C,x0,SSE] = fit_lorentz2(data,min_x,min_y,B0,lambda_min,lambda_max,plot_lower,plot_upper);
        else
            waitfor(errordlg('Wrong number of input arguments!','Error in "fit_lorentz.m"','modal'));
        end
    end
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determination of the start values for the fit
%
%       x0 = min_x
%       T_max = L(inf) = C = max(data(:,2))
%       T_min = L(x0) = C - A = min_y => A = C-T_min = fixed !!!
%       B = FWHM/2 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,B,C,x0,SSE] = fit_lorentz1(data,min_x,min_y,B0)

    T_max = max(data(:,2));
    T_min = min_y;

    x0=min_x;
    C0=T_max;
    
    Clower=0;
    Blower=0;
    
    Cupper=C0*3/2;
    Bupper=5*B0;

    start_values = [B0,C0];
    lower_bound = [Blower,Clower];
    upper_bound = [Bupper,Cupper];
    options = fitoptions('method','NonlinearLeastSquares','Lower',lower_bound,'Upper',upper_bound,'StartPoint',start_values);
    lorentz=fittype(@(B,C,x) C-(C-T_min)*B^2./((x-x0).^2+B^2),'options',options);

    [res,goodness] = fit(data(:,1),data(:,2),lorentz);
    B = abs(res.B);
    C = res.C;
    A=C-T_min;
    L = @(x) C - A*B^2/((x-x0)^2 + B^2);
    SSE = goodness.sse;

    plot(data(:,1),data(:,2));
    hold('on');
    fplot(L,[data(1,1),data((size(data,1)),1)],'r');
    hold ('off');  
end

function [A,B,C,x0,SSE] = fit_lorentz2(data,min_x,min_y,B0,lambda_min,lambda_max,plot_lower,plot_upper)

    i_lambda_min = lambda_inv(data(:,1)',lambda_min);
    i_lambda_max = lambda_inv(data(:,1)',lambda_max);
    i_plot_lower = lambda_inv(data(:,1)',plot_lower);
    i_plot_upper = lambda_inv(data(:,1)',plot_upper);

    T_max = max(data(:,2));
    T_min = min_y;

    x0=min_x;
    C0=T_max;
    
    Clower=0;
    Blower=0;
    
    Cupper=C0*3/2;
    Bupper=5*B0;

    start_values = [B0,C0];
    lower_bound = [Blower,Clower];
    upper_bound = [Bupper,Cupper];
    options = fitoptions('method','NonlinearLeastSquares','Lower',lower_bound,'Upper',upper_bound,'StartPoint',start_values);   
    lorentz=fittype(@(B,C,x) C-(C-T_min).*B^2./((x-x0).^2+B^2),'options',options);

    [res,goodness] = fit(data(i_lambda_min:i_lambda_max,1),data(i_lambda_min:i_lambda_max,2),lorentz);
    B = abs(res.B);
    C = res.C;
    A = C-T_min;
    L = @(x) C - A*B^2/((x-x0)^2 + B^2);
    SSE = goodness.sse;
    
    plot(data(i_plot_lower:i_plot_upper,1),data(i_plot_lower:i_plot_upper,2));
    hold('on');
    fplot(L,[plot_lower,plot_upper],'r');
    hold('off');

end

