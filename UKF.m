fid = fopen('result_mse.txt', 'w');
for number = 0:0
    Q = eye(6);
    R = [1];
    H = [0, 0, 0, 0, 1, 0];
    a = 1.5;

    filename = strrep('data/xlsx/SNR1.xlsx', '1', num2str(number));
    sheet = 1;
    time = xlsread(filename, sheet, 'A:A');   
    output = xlsread(filename, sheet, 'B:B');
    input = xlsread(filename, sheet, 'C:C');
    ideal = xlsread(filename, sheet, 'D:D');
    
    n_samples = 700;
    [row, column] = size(output);
    interval = row / n_samples;
    indices = int32(1:interval:row);
    if length(indices) < n_samples
        indices = [indices row];
    end
    time = time(indices);
    input = input(indices);
    output = output(indices);
    ideal = ideal(indices);

    delta = diff(time);
    input = input(2:end);
    output = output(2:end);

    mu = [0; 0; 0; 0; 0; 0];
    P = 10 + sqrt(1) * randn(6, 6);
    P(P<0) = -P(P<0);

    anse= [];
    [row, column] = size(output);
    for i = 1:row
        X = iterate(mu, P, a);
        [X, Y] = F(X, input(i), delta(i));
        X = real(X);
        Y = real(Y);

        [x_pred, p_pred] = predict(X, a, Q);
        [mu, P] = update(x_pred, p_pred, output(i), H, R);

        anse= [anse Y(1)];
    end

    error = rmse(anse.',ideal(1:699)) * 1.0;
    fprintf(fid, "RMSE for %s is %f\n", filename, error);
    fprintf("RMSE for %s is %f\n", filename, error);

    plot(anse, 'Color', "blue");
    hold on;
    plot(ideal, 'Color', "red");
    hold on;
    plot(input, 'LineWidth', 2, 'Color', "green");
    legend('result', 'ideal_output', 'input');
    filename = strrep(filename, '.xlsx', '.fig');
    filename = strrep(filename, 'xlsx', 'fig');
    savefig(filename);
    clf;
    close;
end

fclose(fid);

function [mu, P, y] = update(x_pred, p_pred, y_real, H, R)
    y = H * x_pred;
    pyy = H * p_pred * H.' + R;
    pxy = p_pred * H.';
    x = x_pred + pxy * inv(pyy) * (y_real - y);
    P = p_pred - pxy * inv(pyy) * pyy.';
    mu = mean(x, 2);
end


function [x_pred, p_pred] = predict(X, a, Q)
    mu = mean(X, 2);
    L = 6;
    w = [(1-1/a^2)];
    for i = 1:2*L
        w = [w; (1/2/L/a/a)];
    end
    x_pred = ((w.') * (X.')).';
    p_pred = zeros(L, L);
    for i = 1:(2*L+1)
        tmp = (X(:, i) - x_pred) * (X(:, i) - x_pred).';
        p_pred = p_pred + w(i) * tmp;
    end
    p_pred = p_pred + Q;
end

function X = iterate(mu, P, a)
    L = 6;
    
    X = [mu];
    sq_root = transpose(sqrtm(L * P));
    %sq_root
    for i = 1:L
        tmp = mu + a*sq_root(:, i);
        X = [X, tmp];
    end
    for i = 1:L
        tmp = mu - a*sq_root(:, i);
        X = [X, tmp];
    end
end

function [Xnew, Y] = F(Xold, u, dT)
    Xnew = [];
    Y = [];
    %Xold
    for X = Xold
        w = 2 * pi * (1/15);
        m1=1; %元件有無noise
        m2=0;
        m3=0;
        m4=0;
        m5=0;
        Zr2=1000000; %1M
        Zr1=1000000; %1M
        Zc1=1i * w* (10^-6);%w=?
        Zc2=1i * w* (10^-6);
        N1 = X(1);
        N2 = X(2);
        N3 = X(3);
        N4 = X(4);
        
        A55= -1/((Zr1+m1*N1)*(Zc1+m3*N3))-1/((Zr2+m2*N2)*(Zc1+m3*N3));
        A56= 1/((Zr2+m2*N2)*(Zc1+m3*N3));
        A65= 1/((Zr2+m2*N2)*(Zc2+m4*N4));
        A66= -1/((Zr2+m2*N2)*(Zc2+m4*N4));
        A= [1-dT 0 0 0 0 0;
            0 1-dT 0 0 0 0;
            0 0 1-dT 0 0 0;
            0 0 0 1-dT 0 0;
            0 0 0 0 1+dT*A55 A56;
            0 0 0 0 A65 1+dT*A66]; %𝜌= 1

        Z5= dT*u/((Zr1+m1*N1) * (Zc1+m3*N3));
        Z= transpose([0 0 0 0 Z5 0]);

        K56= m5/((Zr1+m1*N1) * (Zc1+m3*N3)); %k該是matrix還是一個值？
        K= [0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 K56 0;
            0 0 0 0 0 0]; %非常不確定
        B= transpose([wgn(1, 5, 10) 0]); % power of noise 50dbW, brownian motion process微分是高斯雜訊吧
        c=[0 0 0 0 1 0];

        Xtmp= A*X+ Z+ K*B;
        Ytmp=c*X;
        Xnew = [Xnew Xtmp];
        Y = [Y Ytmp];
    end
end

function error = rmse(a, b)
    error = 0;
    [rowa, columna] = size(a);
    for i = 1:rowa
        x = a(i);
        y = b(i);
        error = error + ((x - y) ^ 2) / rowa;
    end
    error = sqrt(error);
end