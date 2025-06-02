%
%    some test  on Articulation index, HL0dB
%    Irino, T.
%   Created: 8 Feb 22
%   Modified: 8 Feb 22
%
%   AI from 
%     Daniel R. Raichel, "THE SCIENCE AND APPLICATIONS OF ACOUSTICS (2nd Ed.)",
%     Springer. 2006
%

NumOneThird = 7*3+1;
FreqOneThird = exp(linspace(log(125),log(8000),NumOneThird));
NumOneSixth = 7*6+1;
FreqOneSixth = exp(linspace(log(125),log(8000),NumOneSixth));
[~, n1000] = min(abs(FreqOneSixth-1000));
Bw1000 = FreqOneSixth(n1000 + 1) - FreqOneSixth(n1000 - 1);

AIfreq = [200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000];
AIweight = [4,  10 ,  10,   14,   14,   20,   20,    24,     30,     37,    37,     34,     34,    24,     20];

ERBband = linspace(Freq2ERB(100),Freq2ERB(8000),100);
Fr1 = ERB2Freq(ERBband);
[~, ERBw] = Freq2ERB(1000);
AIweightERB = interp1(AIfreq,AIweight,Fr1,'linear','extrap')* ERBw/Bw1000;
pwr = 0.5;
AIweightERB  = max(AIweightERB,0).^pwr*5;
[Fr1; AIweightERB]

SSIparam.Fr1 = Fr1;
SSIparam.F0_limit = 110;
[SSIweight] = F0limit2SSIweight(SSIparam);

subplot(2,1,1)
semilogx(AIfreq,AIweight,'o-', Fr1,AIweightERB,Fr1,SSIweight*10)
grid on
xlabel('Frequency (Hz)')
ylabel('AI weight (ERB interp. pwr0.3)');


subplot(2,1,2)
plot(Freq2ERB(AIfreq),AIweight,'o-', ERBband,AIweightERB,ERBband,SSIweight*10)
grid on
xlabel('ERB_N number (Cam)')
ylabel('AI weight (ERB interp. pwr0.3)');

find(Fr1 < 200)

