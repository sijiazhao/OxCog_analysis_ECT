function y = fconv(x,convWidth,kernel)

% kernel = 'uniform';

switch kernel
    case 'poisson'
        convWidth = 100;
        Fs = 5;
        a = 1/Fs;
        
        convWindow = zeros(convWidth, 1);
        for i = 1:length(convWindow)
            convWindow(i) = (a^2)*i*exp(-a*i);
            %     Window(i) = sqrt(1/a/pi)*exp(-(1/a)*((i-3*(1/a))/Fs)^2);
            %     Window(i) = 1/2*(1-cos(2*pi*i/(Fs-1)))*exp(-a*abs(Fs-1-2*i)/(Fs-1));
        end
        convWindow = convWindow/sum(convWindow);
        %         plot(convWindow)
        y = conv2(x,convWindow); %lucas code
        % y = filter(convWindow,1,x); %sijia code
        
    case 'uniform' %rectangular window
        
        convWindow = 1/convWidth*ones(convWidth,1);
        y = conv2(x,convWindow); %lucas code
        
    case 'gaussian'
        t = 1:1:convWidth;
        convWindow = gaussmf(t,[convWidth/5 convWidth/2]);
        convWindow = convWindow/sum(convWindow);
        convWindow = convWindow';
        %         plot(t,convWindow)
        y = conv2(x,convWindow); %lucas code
end
end

%% Plot an example of convolution's result
% figure(1000); clf;
% subplot(2,1,1);
% plot(convWindow,'k');
% xlabel('[ms]');
% title(['convolution filter with convWidth=', num2str(convWidth)]);
%
% subplot(2,1,2);
% hold on;
% plot(x,'--k');
% plot(y,'r');
% hold off;
%
% xlabel('time from onset [ms]');
% legend({'original';'convoluted'});
% title(['convolution applied before & after']);
%
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8]);
% filename = ['eg_conv'];
% saveas(gcf,['./Plots/' filename, '.png'],'png');
