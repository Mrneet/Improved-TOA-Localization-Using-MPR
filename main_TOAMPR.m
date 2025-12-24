% TOA-MPR, for SNR and range
% s1 is the reference
clear all;close all;clc;

% Simulation Configurations ***************************************
noiseOrRange = 'nse';
% noiseOrRange = 'rag';
% *****************************************************************


senPos = [ 0       0       0
    3.46    4.44    1.55
    3.05   -2.94    1.84
    0.99    2.92   -0.84
    -0.74    4.80   -0.35
    0.45   -0.09    2.17
    -4.72   -0.73    2.91]';
theta = 22.13*pi/180;
phi = 14.41*pi/180;
u0 = [cos(theta)*cos(phi); sin(theta)*cos(phi); sin(phi)];

[N,M] = size(senPos);

switch noiseOrRange
    case 'nse'
        %******* vs. noise power config *******
        nsePwr = -100:10:40;  % 10log(*)
        % souRange = 80;        % m, far-field
        souRange = 5;         % m, near-field
    case 'rag'
        %******* vs. range config *******
        % nsePwr = -10;  % 10log(*)
        % range = [5,10,20:20:200];  % m
        % range = [10,40:40:320];
end

mon = 500;

rng('Default');
nseTmp = zeros(M,mon);
for m = 1:mon
    nseTmp(:,m) = randn(M,1);
end
nseTmp = nseTmp - mean(nseTmp,2);
nAlg = 6;

totalTime = zeros(nAlg,1);
OPTIONS = optimset('Algorithm','Levenberg-Marquardt','Display','off');

for ir = 1:length(souRange)
    g = 1/souRange(ir);
    souLoc = souRange(ir)*u0;
    r = sqrt(sum((souLoc-senPos).^2,1))';

    for in = 1:length(nsePwr)
        %         disp(['noise: ', num2str(nsePwr(in)),' ...']);
        disp(['range: ',num2str(souRange(ir)),', noise: ', num2str(nsePwr(in)),' ...']);
        Q = 10^(nsePwr(in)/10) * eye(M);
        Qt = 10^(nsePwr(in)/10) * 0.5*(eye(M-1)+ones(M-1));

        % CRLB
        [CRB_MPR,CRB_Cart] = CRLB_TOA(senPos,souLoc,Q);
        crlb_t(ir,in) = CRB_MPR(1,1);
        crlb_p(ir,in) = CRB_MPR(2,2);
        crlb_a(ir,in) = trace(CRB_MPR(1:2,1:2));
        crlb_g(ir,in) = CRB_MPR(3,3);

        crlb_u(ir,in) = trace(CRB_Cart);

        mprSol = zeros(N,mon,nAlg);
        carSol = zeros(N,mon,nAlg);
        for m = 1:mon
            nse = sqrtm(Q)*nseTmp(:,m);
            rNsed = r + nse;

            ia = 1;
            %%%% proposed 2S-SDR
            [mpr2SDR,car2SDR] = TOA_MPR_2SSDR(senPos,rNsed,Q);
            totalTime(ia) = totalTime(ia) + toc;
            mprSol(:,m,ia) = mpr2SDR;
            carSol(:,m,ia) = car2SDR;

            %%%% Newton-Guass iteration
            ia = ia + 1; tic;
            [mprMLE,carMLE] = TOA_MPR_MLE(iniMPR,senPos,rNsed,Q);
            totalTime(ia) = totalTime(ia) + toc;
            mprSol(:,m,ia) = mprMLE;
            carSol(:,m,ia) = carMLE;

            %%%% 2WLS in Cartesian
            ia = ia + 1; tic;
            [x_ls,mpr2WLS] =TOA_2WLS_Cart(senPos,rNsed,Q);
            totalTime(ia) = totalTime(ia) + toc;
            carSol(:,m,ia) = x_ls;
            mprSol(:,m,ia) = mpr2WLS;
            

            %%%% MLE in Cartesian
            ia = ia + 1; tic;
            x_ml0 = lsqnonlin(@(x)fun_toa_MLE(x,rNsed,senPos,Q),x_ls,[],[], OPTIONS );
            x_ml = sign(real(x_ml0)).*abs(x_ml0);
            totalTime(ia) = totalTime(ia) + toc;
            carSol(:,m,ia) = x_ml;
            mprSol(:,m,ia) = [atan2(x_ml(2),x_ml(1)); atan2(x_ml(3),norm(x_ml(1:2))); 1/norm(x_ml)];
        end
        nAlg = ia;
        
        for ia = 1:nAlg
            mse_t(ir,in,ia) = mean((mprSol(1,:,ia)-theta).^2,2);
            mse_p(ir,in,ia) = mean((mprSol(2,:,ia)-phi).^2,2);
            mse_a(ir,in,ia) = mean((mprSol(1,:,ia)-theta).^2,2) + mean((mprSol(2,:,ia)-phi).^2,2);
            mse_g(ir,in,ia) = mean((mprSol(3,:,ia)-g).^2,2);
            mse_u(ir,in,ia) = mean(sum((carSol(:,:,ia)-souLoc).^2,1));
            bias_t(ir,in,ia) = mean(mprSol(1,:,ia),2)-theta;
            bias_p(ir,in,ia) = mean(mprSol(2,:,ia),2)-phi;
            bias_a(ir,in,ia) = sqrt(abs(bias_t(ir,in,ia)).^2 + abs(bias_p(ir,in,ia)).^2);
            bias_g(ir,in,ia) = abs(mean(mprSol(3,:,ia),2)-g);
            % bias_u(ir,in,ia) = sum(abs(mean(carSol(:,:,ia)-souLoc,2)));
            bias_u(ir,in,ia) = mean(sqrt(sum((carSol(:,:,ia)-souLoc).^2,1)));
        end
    end
end

names = {'SDR','2S-SDR','MLE','WLS-Cart','MLE-Cart'};
symbs = 'ox^s+dvp<>';
switch noiseOrRange
    case 'nse'
        xlabtext = '10log(\sigma_r^2(m^2))';
        xdata = nsePwr;
        fileName = ['mat_3D_TOAMPR_Noise_Range',num2str(souRange),'_',datestr(now,'yyyymmmdd_HHMM'),'.mat'];
        % yl_mse = [-90, 15;
        %     -135,0];
        % yl_bias = [-185,10;
        %     -230,-0];
        mstr = 'Noise';
    case 'rag'
        xlabtext = 'Range(m)';
        xdata = souRange;
        fileName = ['mat_3D_TOAMPR_Range_Noise',num2str(nsePwr),'_',datestr(now,'yyyymmmdd_HHMM'),'.mat'];
        % yl_mse = [-40,10;
        %     -83,-40];
        % yl_bias = [-85,0;
        %     -150,-45];
        mstr = 'range';
end
save(fileName);
system('shutdown -s -t 60')
    
% MSE
figure;
subplot(2,1,1)
for ii = 1:nAlg
    plot(xdata, 10*log10(reshape(mse_a(:,:,ii),[],1)), ...
        symbs(ii), 'LineWidth', 1.5, 'DisplayName', names{ii});hold on;grid on;
end
plot(xdata, 10*log10(crlb_a), '--', 'LineWidth', 1.5, 'DisplayName', 'CRLB');
xlabel(xlabtext, 'FontSize', 11);
ylabel('10log(MSE(\theta,\phi)(rad^2))', 'FontSize', 11);
legend('Show', 'FontSize',11, 'Location', 'Northwest', 'NumColumns', 2);

subplot(2,1,2);
for ii = 1:nAlg
    plot(xdata, 10*log10(reshape(mse_g(:,:,ii),[],1)), ...
        symbs(ii), 'LineWidth', 1.5, 'DisplayName', names{ii});hold on;grid on;
end
plot(xdata, 10*log10(crlb_g), '--', 'LineWidth', 1.5, 'DisplayName', 'CRLB');
xlabel(xlabtext, 'FontSize', 11);
ylabel('10log(MSE(g)(m^{-2}))', 'FontSize', 11);
legend('Show', 'FontSize',11, 'Location', 'Northwest', 'NumColumns', 2);
set(gcf,'Position', [680 380 560 500]);

% bias
figure;
subplot(2,1,1)
for ii = 1:nAlg
    plot(xdata, 20*log10(reshape(bias_a(:,:,ii),[],1)), ...
        symbs(ii), 'LineWidth', 1.5, 'DisplayName', names{ii});hold on;grid on;
end
% plot(xdata, 10*log10(crlb_a), '-', 'LineWidth', 1.5, 'DisplayName', 'CRLB-MPR');
xlabel(xlabtext, 'FontSize', 11);
ylabel('10log(Bias(\theta,\phi)(rad^2))', 'FontSize', 11);
legend('Show', 'FontSize',11, 'Location', 'Northwest', 'NumColumns', 2);

subplot(2,1,2);
for ii = 1:nAlg
    plot(xdata, 20*log10(reshape(bias_g(:,:,ii),[],1)), ...
        symbs(ii), 'LineWidth', 1.5, 'DisplayName', names{ii});hold on;grid on;
end
% % plot(xdata, 10*log10(crlb_g), '-', 'LineWidth', 1.5, 'DisplayName', 'CRLB-MPR');
xlabel(xlabtext, 'FontSize', 11);
ylabel('10log(Bias(g)(m^{-2}))', 'FontSize', 11);
legend('Show', 'FontSize',11, 'Location', 'Northwest', 'NumColumns', 2);
set(gcf,'Position', [680 380 560 500]);

% (MSE of MPR) / CRLB
figure;
subplot(2,1,1);
for ii = 1:nAlg
    plot(xdata,10*log10(reshape(mse_a(:,:,ii),[],1)./reshape(crlb_a,[],1)),symbs(ii),'linewidth',1.5,'DisplayName',names{ii});hold on;grid on;
end
xlabel('10log(\sigma_n^2(m^2))','FontSize',11);ylabel('10log(MSE(\theta,\phi)/CRLB(\theta,\phi))','FontSize',11);
legend('show','Location','Northwest','FontSize',11,'NumColumns',3);
ylim([-2 5])

subplot(2,1,2);
for ii = 1:nAlg
    plot(xdata, 10*log10(reshape(mse_g(:,:,ii),[],1)./reshape(crlb_g,[],1)), ...
        symbs(ii), 'LineWidth', 1.5, 'DisplayName', names{ii});hold on;grid on;
end
xlabel('10log(\sigma_n^2(m^2))','FontSize',11);ylabel('10log(MSE(g)/CRLB(g))','FontSize',11);
legend('show','Location','Northwest','FontSize',11,'NumColumns',3);
ylim([-2 5])
set(gcf,'Position', [680 380 560 500]);


% if length(nsePwr)>1
%     %     load('mat_3D_TOAMPR_Noise_Range80_2018Dec10_2316.mat','mse_t1','mse_t2','mse_t3',...
%     %     'mse_p1','mse_p2','mse_p3','mse_g1','mse_g2','mse_g3','mse_u1','mse_u2','mse_u3');
%     %     fileName = ['mat_3D_TOAMPR_Noise_',datestr(now,'yyyymmmdd_HHMM'),'.mat'];
%     fileName = ['mat_3D_TOAMPR_Noise_Range',num2str(range),'_',datestr(now,'yyyymmmdd_HHMM'),'.mat'];
%     save(fileName);
% 
%     figure;
%     subplot(3,1,1);
%     plot(nsePwr,10*log10(mse_t1),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     plot(nsePwr,10*log10(mse_t2),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     plot(nsePwr,10*log10(mse_t3),'^','linewidth',1.5,'DisplayName','MLE');
%     plot(nsePwr,10*log10(mse_t4),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     plot(nsePwr,10*log10(mse_t5),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     %     plot(nsePwr,10*log10(mse_t6),'d','linewidth',1.5,'DisplayName','Algebric');
%     plot(nsePwr,10*log10(crlb_t),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('10log(\sigma_r^2(m^2))','FontSize',11);ylabel('10log(MSE(\theta)(rad^2))','FontSize',11);
%     leg1 = legend('show','Location','Northwest');
%     set(leg1,'FontSize',11);
% 
%     subplot(3,1,2);
%     plot(nsePwr,10*log10(mse_p1),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     plot(nsePwr,10*log10(mse_p2),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     plot(nsePwr,10*log10(mse_p3),'^','linewidth',1.5,'DisplayName','MLE');
%     plot(nsePwr,10*log10(mse_p4),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     plot(nsePwr,10*log10(mse_p5),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     %     plot(nsePwr,10*log10(mse_p6),'d','linewidth',1.5,'DisplayName','Algebric');
%     plot(nsePwr,10*log10(crlb_p),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('10log(\sigma_r^2(m^2))','FontSize',11);ylabel('10log(MSE(\phi)(rad^2))','FontSize',11);
%     leg1 = legend('show','Location','Northwest');
%     set(leg1,'FontSize',11);
% 
%     subplot(3,1,3);
%     plot(nsePwr,10*log10(mse_g1),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     plot(nsePwr,10*log10(mse_g2),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     plot(nsePwr,10*log10(mse_g3),'^','linewidth',1.5,'DisplayName','MLE');
%     plot(nsePwr,10*log10(mse_g4),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     plot(nsePwr,10*log10(mse_g5),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     %     plot(nsePwr,10*log10(mse_g6),'d','linewidth',1.5,'DisplayName','Algebric');
%     plot(nsePwr,10*log10(crlb_g),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('10log(\sigma_r^2(m^2))','FontSize',11);ylabel('10log(MSE(g)(m^{-2}))','FontSize',11);
%     leg1 = legend('show','Location','Northwest');
%     set(leg1,'FontSize',11);
%     set(gcf,'Position', [680 180 560 750])
% 
%     % MSE of angle and inverse-range
%     figure;
%     subplot(2,1,1);
%     plot(nsePwr,10*log10(mse_t1+mse_p1),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     plot(nsePwr,10*log10(mse_t2+mse_p2),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     plot(nsePwr,10*log10(mse_t3+mse_p3),'^','linewidth',1.5,'DisplayName','MLE');
%     plot(nsePwr,10*log10(mse_t4+mse_p4),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     plot(nsePwr,10*log10(mse_t5+mse_p5),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     plot(nsePwr,10*log10(crlb_t+crlb_p),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('10log(\sigma_n^2(m^2))','FontSize',11);ylabel('10log(MSE(\theta,\phi)(rad^2))','FontSize',11);
%     leg1 = legend('show','Location','Northwest');
%     set(leg1,'FontSize',11);
%     ylim([-120 40])
% 
%     subplot(2,1,2);
%     plot(nsePwr,10*log10(mse_g1),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     plot(nsePwr,10*log10(mse_g2),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     plot(nsePwr,10*log10(mse_g3),'^','linewidth',1.5,'DisplayName','MLE');
%     plot(nsePwr,10*log10(mse_g4),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     plot(nsePwr,10*log10(mse_g5),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     plot(nsePwr,10*log10(crlb_g),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('10log(\sigma_n^2(m^2))','FontSize',11);ylabel('10log(MSE(g)(m^{-2}))','FontSize',11);
%     leg1 = legend('show','Location','Northwest');
%     set(leg1,'FontSize',11);
%     set(gcf,'Position', [680 380 560 580]);
% 
%     % (MSE of MPR) / CRLB
%     figure;
%     subplot(2,1,1);
%     plot(nsePwr,10*log10((mse_t1+mse_p1)./(crlb_t+crlb_p)),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     plot(nsePwr,10*log10((mse_t2+mse_p2)./(crlb_t+crlb_p)),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     plot(nsePwr,10*log10((mse_t3+mse_p3)./(crlb_t+crlb_p)),'^','linewidth',1.5,'DisplayName','MLE');
%     plot(nsePwr,10*log10((mse_t4+mse_p4)./(crlb_t+crlb_p)),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     plot(nsePwr,10*log10((mse_t5+mse_p5)./(crlb_t+crlb_p)),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     %     plot(nsePwr,10*log10(crlb_t+crlb_p),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('10log(\sigma_n^2(m^2))','FontSize',11);ylabel('10log(MSE(\theta,\phi)/CRLB(\theta,\phi))','FontSize',11);
%     leg1 = legend('show','Location','Northwest');
%     set(leg1,'FontSize',11);
%     ylim([-2 10])
% 
%     subplot(2,1,2);
%     plot(nsePwr,10*log10(mse_g1./crlb_g),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     plot(nsePwr,10*log10(mse_g2./crlb_g),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     plot(nsePwr,10*log10(mse_g3./crlb_g),'^','linewidth',1.5,'DisplayName','MLE');
%     plot(nsePwr,10*log10(mse_g4./crlb_g),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     %     plot(nsePwr,10*log10(mse_g5./crlb_g),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     %     plot(nsePwr,10*log10(crlb_g),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('10log(\sigma_n^2(m^2))','FontSize',11);ylabel('10log(MSE(g)/CRLB(g))','FontSize',11);
%     leg1 = legend('show','Location','Northwest');
%     set(leg1,'FontSize',11);
%     %     set(gcf,'Position', [680 380 560 580]);
% 
%     %     export_fig -transparent mse_nseLvl_range80.pdf
% 
%     % (MSE of u) / CRLB
%     figure;
%     plot(nsePwr,10*log10(mse_u1./crlb_u),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     plot(nsePwr,10*log10(mse_u2./crlb_u),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     plot(nsePwr,10*log10(mse_u3./crlb_u),'^','linewidth',1.5,'DisplayName','MLE');
%     plot(nsePwr,10*log10(mse_u4./crlb_u),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     plot(nsePwr,10*log10(mse_u5./crlb_u),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     %     plot(nsePwr,10*log10(mse_u6./crlb_u),'d','linewidth',1.5,'DisplayName','Algebric');
%     %     plot(nsePwr,10*log10(crlb_u),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('10log(\sigma_n^2(m^2))','FontSize',11);ylabel('10log(MSE(u)/CRLB)','FontSize',11);
%     leg1 = legend('show','Location','Northwest');
%     set(leg1,'FontSize',11);
%     ylim([-2 10])
%     set(gcf,'Position', [680 380 560 300]);
% 
%     %     export_fig -transparent mse_nseLvl_posRatio_r80.pdf
% elseif length(range)>1
%     %     load('mat_3D_TOAMPR_Range_Noise-20_2018Dec12_1513.mat','mse_t1','mse_t2','mse_t3',...
%     %     'mse_p1','mse_p2','mse_p3','mse_g1','mse_g2','mse_g3','mse_u1','mse_u2','mse_u3');
%     fileName = ['mat_3D_TOAMPR_Range_Noise',num2str(nsePwr),'_',datestr(now,'yyyymmmdd_HHMM'),'.mat'];
%     save(fileName);
% 
%     figure;
%     subplot(3,1,1);
%     plot(range,10*log10(mse_t1),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     plot(range,10*log10(mse_t2),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     plot(range,10*log10(mse_t3),'^','linewidth',1.5,'DisplayName','MLE');
%     plot(range,10*log10(mse_t4),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     plot(range,10*log10(mse_t5),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     %     plot(range,10*log10(mse_t6),'d','linewidth',1.5,'DisplayName','Algebric');
%     plot(range,10*log10(crlb_t),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('10log(\sigma_r^2(m^2))','FontSize',11);ylabel('10log(MSE(\theta)(rad^2))','FontSize',11);
%     leg1 = legend('show','Location','Northwest');
%     set(leg1,'FontSize',11);
% 
%     subplot(3,1,2);
%     plot(range,10*log10(mse_p1),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     plot(range,10*log10(mse_p2),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     plot(range,10*log10(mse_p3),'^','linewidth',1.5,'DisplayName','MLE');
%     plot(range,10*log10(mse_p4),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     plot(range,10*log10(mse_p5),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     %     plot(range,10*log10(mse_p6),'d','linewidth',1.5,'DisplayName','Algebric');
%     plot(range,10*log10(crlb_p),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('10log(\sigma_r^2(m^2))','FontSize',11);ylabel('10log(MSE(\phi)(rad^2))','FontSize',11);
%     leg1 = legend('show','Location','Northwest');
%     set(leg1,'FontSize',11);
% 
%     subplot(3,1,3);
%     plot(range,10*log10(mse_g1),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     plot(range,10*log10(mse_g2),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     plot(range,10*log10(mse_g3),'^','linewidth',1.5,'DisplayName','MLE');
%     plot(range,10*log10(mse_g4),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     plot(range,10*log10(mse_g5),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     %     plot(nsePwr,10*log10(mse_g6),'d','linewidth',1.5,'DisplayName','Algebric');
%     plot(range,10*log10(crlb_g),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('10log(\sigma_r^2(m^2))','FontSize',11);ylabel('10log(MSE(g)(m^{-2}))','FontSize',11);
%     leg1 = legend('show','Location','Northwest');
%     set(leg1,'FontSize',11);
%     set(gcf,'Position', [680 180 560 750])
% 
%     figure;
%     subplot(2,1,1);
%     plot(range,10*log10(mse_t1+mse_p1),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     plot(range,10*log10(mse_t2+mse_p2),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     plot(range,10*log10(mse_t3+mse_p3),'^','linewidth',1.5,'DisplayName','MLE');
%     plot(range,10*log10(mse_t4+mse_p4),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     plot(range,10*log10(mse_t5+mse_p5),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     plot(range,10*log10(crlb_t+crlb_p),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('Range (m)','FontSize',11);ylabel('10log(MSE(\theta,\phi)(rad^2))','FontSize',11);
%     leg1 = legend('show','Location','Northeast');
%     set(leg1,'FontSize',11);
% 
%     subplot(2,1,2);
%     plot(range,10*log10(mse_g1),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     plot(range,10*log10(mse_g2),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     plot(range,10*log10(mse_g3),'^','linewidth',1.5,'DisplayName','MLE');
%     plot(range,10*log10(mse_g4),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     plot(range,10*log10(mse_g5),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     plot(range,10*log10(crlb_g),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('Range (m)','FontSize',11);ylabel('10log(MSE(g)(m^{-2}))','FontSize',11);
%     leg1 = legend('show','Location','Northeast');
%     set(leg1,'FontSize',11);
%     set(gcf,'Position', [680 380 560 500]);
% 
%     %     export_fig -transparent mse_range_noise-20.pdf
% 
%     figure;
%     plot(range,10*log10(mse_u1),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     plot(range,10*log10(mse_u2),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     plot(range,10*log10(mse_u3),'^','linewidth',1.5,'DisplayName','MLE');
%     plot(range,10*log10(mse_u4),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     plot(range,10*log10(mse_u5),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     plot(range,10*log10(crlb_u),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('Range (m)','FontSize',11);ylabel('10log(MSE(u)(m^{2}))','FontSize',11);
%     leg1 = legend('show','Location','Southeast');
%     set(leg1,'FontSize',11);
%     set(gcf,'Position', [680 380 560 300]);
% 
%     figure;
%     plot(range,10*log10(mse_u1./crlb_u),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     plot(range,10*log10(mse_u2./crlb_u),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     plot(range,10*log10(mse_u3./crlb_u),'^','linewidth',1.5,'DisplayName','MLE');
%     plot(range,10*log10(mse_u4./crlb_u),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     plot(range,10*log10(mse_u5./crlb_u),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     xlabel('Range (m)','FontSize',11);ylabel('10log(MSE(u)/CRLB)','FontSize',11);
%     leg1 = legend('show','Location','Northwest');
%     set(leg1,'FontSize',11);
% 
%     figure;
%     subplot(3,1,1);
%     p(1) = plot(range,10*log10(mse_t1+mse_p1),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     p(2) = plot(range,10*log10(mse_t2+mse_p2),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     p(3) = plot(range,10*log10(mse_t3+mse_p3),'^','linewidth',1.5,'DisplayName','MLE');
%     p(4) = plot(range,10*log10(mse_t4+mse_p4),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     p(5) = plot(range,10*log10(mse_t5+mse_p5),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     p(6) = plot(range,10*log10(crlb_t+crlb_p),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('Range (m)','FontSize',11);ylabel('10log(MSE(\theta,\phi)(rad^2))','FontSize',11);
%     %     leg1 = legend('show','Location','Northeast','Orientation','horizontal');
%     ah=axes('position',get(gca,'position')+[-0.2 0 0 0],'visible','off');
%     leg1 = legend(ah,p(1:3),'SDR','2S-SDR','MLE');
%     set(leg1,'FontSize',10);
%     ah=axes('position',get(gca,'position')+[+0.2 0 0 0],'visible','off');
%     leg1 = legend(ah,p(4:6),'WLS-Cart','MLE-Cart','CRLB');
%     set(leg1,'FontSize',10);
% 
%     subplot(3,1,2);
%     p(1) = plot(range,10*log10(mse_g1),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     p(2) = plot(range,10*log10(mse_g2),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     p(3) = plot(range,10*log10(mse_g3),'^','linewidth',1.5,'DisplayName','MLE');
%     p(4) = plot(range,10*log10(mse_g4),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     p(5) = plot(range,10*log10(mse_g5),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     p(6) = plot(range,10*log10(crlb_g),'--','linewidth',1.5,'DisplayName','CRLB');
%     xlabel('Range (m)','FontSize',11);ylabel('10log(MSE(g)(m^{-2}))','FontSize',11);
%     %     leg1 = legend('show','Location','Northeast');
%     %     set(leg1,'FontSize',10);
%     ah=axes('position',get(gca,'position')+[-0.2 0 0 0],'visible','off');
%     leg1 = legend(ah,p(1:3),'SDR','2S-SDR','MLE');
%     set(leg1,'FontSize',10);
%     ah=axes('position',get(gca,'position')+[+0.2 0 0 0],'visible','off');
%     leg1 = legend(ah,p(4:6),'WLS-Cart','MLE-Cart','CRLB');
%     set(leg1,'FontSize',10);
% 
%     subplot(3,1,3)
%     p(1) = plot(range,10*log10(mse_u1./crlb_u),'o','linewidth',1.5,'DisplayName','SDR');hold on;grid on;
%     p(2) = plot(range,10*log10(mse_u2./crlb_u),'x','linewidth',1.5,'DisplayName','2S-SDR');
%     p(3) = plot(range,10*log10(mse_u3./crlb_u),'^','linewidth',1.5,'DisplayName','MLE');
%     p(4) = plot(range,10*log10(mse_u4./crlb_u),'s','linewidth',1.5,'DisplayName','WLS-Cart');
%     p(5) = plot(range,10*log10(mse_u5./crlb_u),'+','linewidth',1.5,'DisplayName','MLE-Cart');
%     xlabel('Range (m)','FontSize',11);ylabel('10log(MSE(u)/CRLB(u))','FontSize',11);
%     %     leg1 = legend('show','Location','Northeast');
%     %     set(leg1,'FontSize',10);
%     leg1 = legend(p(1:3),'SDR','2S-SDR','MLE','Location','Northwest');
%     set(leg1,'FontSize',10);
%     ah=axes('position',get(gca,'position')+[0.18 0 0 0],'visible','off');
%     leg1 = legend(ah,p(4:5),'WLS-Cart','MLE-Cart','Location','Northwest');
%     set(leg1,'FontSize',10);
%     set(gcf,'Position', [680 324 560 650]);
% 
%     %     export_fig -transparent mse&ratio_range_noise-10.pdf
% end


function output = fun_toa_MLE(x,r,s,Q)
    output = (r-sqrt(sum((x-s).^2,1))')'/Q*(r-sqrt(sum((x-s).^2,1))');
end
