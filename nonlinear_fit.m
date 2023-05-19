clc;
clear;
close all;
exp_LTP_raw=[   0	2.35E-11
1	2.35E-11
2	9.62E-11
3	1.93E-10
4	3.08E-10
5	4.98E-10
6	6.64E-10
7	8.48E-10
8	9.98E-10
9	1.14E-09
10	1.32E-09
11	1.53E-09
12	1.72E-09
13	1.70E-09
14	1.93E-09
15	1.98E-09
16	2.18E-09
17	2.35E-09
18	2.41E-09
19	2.50E-09
20	2.63E-09
21	2.66E-09
22	2.87E-09
23	2.85E-09
24	3.09E-09
25	3.18E-09
26	3.23E-09
27	3.42E-09
28	3.65E-09
29	3.59E-09
30	3.69E-09
31	3.78E-09
32	3.94E-09
33	3.98E-09
34	4.06E-09
35	4.05E-09
36	4.12E-09
37	4.19E-09
38	4.25E-09
39	4.44E-09
40	4.62E-09
41	4.80E-09
42	4.83E-09
43	5.05E-09
44	4.97E-09
45	5.10E-09
46	5.34E-09
47	5.23E-09
48	5.51E-09
49	5.59E-09];

exp_LTD_raw=[   0	1.09E-10
1	1.22E-10
2	1.23E-10
3	1.35E-10
4	1.11E-10
5	1.39E-10
6	1.13E-10
7	1.12E-10
8	1.23E-10
9	1.36E-10
10	1.52E-10
11	1.58E-10
12	1.53E-10
13	1.98E-10
14	1.98E-10
15	2.22E-10
16	2.30E-10
17	2.51E-10
18	2.70E-10
19	2.86E-10
20	3.03E-10
21	3.21E-10
22	3.20E-10
23	3.61E-10
24	4.00E-10
25	4.28E-10
26	4.19E-10
27	4.41E-10
28	4.68E-10
29	5.07E-10
30	5.49E-10
31	6.17E-10
32	6.74E-10
33	6.59E-10
34	7.29E-10
35	8.14E-10
36	9.37E-10
37	1.02E-09
38	1.29E-09
39	1.50E-09
40	1.66E-09
41	1.86E-09
42	2.20E-09
43	2.55E-09
44	2.86E-09
45	3.10E-09
46	3.51E-09
47	3.78E-09
48	4.38E-09
49	5.59E-09];
xf_ltp = floor(max(exp_LTP_raw(:,1)));
xf_ltd = floor(max(exp_LTD_raw(:,1)));
x_step_ltp = 1/xf_ltp;
x_step_ltd = 1/xf_ltd;
yf_ltp = max(exp_LTP_raw(:,2));
yf_ltd = max(exp_LTD_raw(:,2));
yi_ltp = min(exp_LTP_raw(:,2));
yi_ltd = min(exp_LTD_raw(:,2));

exp_LTP(:,1) = exp_LTP_raw(:,1)/xf_ltp;
exp_LTD(:,1) = exp_LTD_raw(:,1)/xf_ltd;
exp_LTP(:,2) = (exp_LTP_raw(:,2)-yi_ltp)/(yf_ltp-yi_ltp);
exp_LTD(:,2) = (exp_LTD_raw(:,2)-yi_ltd)/(yf_ltd-yi_ltd);

plot(exp_LTP(:,1), exp_LTP(:,2), 'bo', 'LineWidth', 1);
hold on;
plot(exp_LTD(:,1), exp_LTD(:,2), 'ro', 'LineWidth', 1);

xf = 1;
A_LTP = 2.5;
B_LTP = 1./(1-exp(-1./A_LTP));
A_LTD = -0.15;
B_LTD = 1./(1-exp(-1./A_LTD));

% LTP fitting
var_amp = 0;    % LTP cycle-to-cycle variation
rng(103);
x_ltp(1) = 0;
y_ltp(1) = 0;
for n=1:1/x_step_ltp+1
    x_ltp(n+1) = x_ltp(n)+x_step_ltp;
    y_ltp(n+1) = B_LTP(1)*(1-exp(-x_ltp(n+1)/A_LTP(1)));
    delta_y = (y_ltp(n+1)-y_ltp(n)) + randn*var_amp;
    y_ltp(n+1) = y_ltp(n) + delta_y;   
    if y_ltp(n+1)>=1
        y_ltp(n+1)=1;
    elseif y_ltp(n+1)<=0
        y_ltp(n+1)=0;
    end
    x_ltp(n+1) = -A_LTP(1)*log(1-(y_ltp(n+1))/B_LTP(1));
end
plot((0:n-1)/(n-1), y_ltp(1:n), 'b', 'linewidth', 2);

% LTD fitting
var_amp = 0;    % LTD cycle-to-cycle variation
rng(898);
x_ltd(1) = 1;
y_ltd(1) = 1;
for n=1:1/x_step_ltd+1
    x_ltd(n+1) = x_ltd(n)-x_step_ltd;
    y_ltd(n+1) = B_LTD(1)*(1-exp(-x_ltd(n+1)/A_LTD(1)));
    delta_y = (y_ltd(n+1)-y_ltd(n)) + randn*var_amp;
    y_ltd(n+1) = y_ltd(n) + delta_y;
    if y_ltd(n+1)>=1
        y_ltd(n+1)=1;
    elseif y_ltd(n+1)<=0
        y_ltd(n+1)=0;
    end
    x_ltd(n+1) = -A_LTD(1)*log(1-(y_ltd(n+1))/B_LTD(1));
end
x_start = numel(x_ltd(:));
x_end = numel(x_ltd(:)) - n;
plot((n-1:-1:0)/(n-1), y_ltd(1:n), 'r', 'linewidth', 2);

xlabel('Normalized Pulse #');
ylabel('Normalized Conductance');
legend('Exp data (LTP)','Exp data (LTD)', 'Fit (LTP)', 'Fit (LTD)', 'location', 'southeast');
