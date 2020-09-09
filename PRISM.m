% Polymer reference interaction site model (PRISM), for modelling of 
% worm-like micelle/polymer solutions. Written by Ashley Williams, 2019-2020.
% 
% Adapted from SASView C-source code.
%

clc
clf
clear
%
% Currently set up with the chi^2 to divide by error:
% chi^2 = Sum[(meas_i-fit_i)^2/err_i^2]
%
Ang = char(197);
%--------------------------------load--------------------------------------
filename = '';%Insert file name (.mat)
load(filename);
data = ;%filename (no quotations)
%---------------------------UI-Parameters----------------------------------

A = 1;            % Overall Scale.
r_cs = 20;        % Cross-sectional Radius
BD = 0.001;       % Background.
L = 1000;         % Contour Length
B = 0.001;          % Peak scale, 0 = no structure factor: NOTE: EITHER MAKE REALLY SMALL OR USE FC MODEL for no Structure factor.
b = 100;          % Kuhn Length 
SLD = 1;          % Scattering length density - Cylinder
SLD_solv = 6.3;   % Scattering length density - solvent
sigma = 1;       % Unknown fit parameter for PRISM
r_eff = 20;       %Effective radius
PD = 0.16;         % Polydispersity, 0 = no polydispersity

%----------------------Schulz-Polydispersity-------------------------------
    qmin = 0.0001;
    qmax = 1;
    qvals = data(:,1)'; % Range of values
    pvals = qvals; 


Z = PD^(-2);
sz = @(x) ((Z+1)^(Z+1))*((x/r_cs)^(Z))*exp(-(Z+1)*(x/r_cs))...
    /(r_cs*gamma(Z+1)); %Schulz distribution function

rvals = (1:1:100); %Number of points for PD calculation.
for j = 1:size(rvals,2)
    dist_r(1,j) = sz(rvals(1,j)); 
end
dist_r = dist_r/max(dist_r); %Normalise the distribution weights.


%---------------------------Other-Constants-------------------------------
contrast = (SLD-SLD_solv)^2; % Del ro^2
volume = 3.1416*r_cs^2*L;

n_b = L./b;

if L > 10*b
    c = @(x) 3.06*(x)^-0.44;
else
    c = @(x) 1;
end

alpha = @(x) sqrt((1+(x/3.12)^2 + (x/8.67)^3)^(0.176/3));

% "Rgsquareshort(L,b),<Rg^2>_0"
% r_g0 = ((L*b)/6)*(1-(3/(2*n_b))+(3/(2*(n_b)^2))-(3/(4*(n_b)^3))*(1-exp(-2*n_b)));

% "Rgsquare(L,b), <Rg^2>"
r_g = (alpha(n_b)^2)*((b*L)/6);

% Rg = <Rg^2>^0.5
r_gsq = sqrt(r_g);
%r_gsq = sqrt((alpha(n_b)^2)*r_g0);

%-------------------------Scattering-functions-----------------------------
%"w_WR(x)
w = @(x) (1+tanh((x-1.523)/0.1477))/2;

u = @(q) r_g*(q^2);

%Debye scattering
S_debye = @(q) 2*(exp(-u(q))+u(q)-1)/u(q)^2;

%Excluded volume scattering
S_exv = @(q) ((1-w(q*r_gsq))*S_debye(q))+(w(q*r_gsq)*(1.22*(q*r_gsq)^(-1/0.585)+0.4288*(q*r_gsq)^(-2/0.585)-1.651*(q*r_gsq)^(-3/0.585)));

%Altered scattering for fixing kink.
S_exv_new = @(q) ((1-w(q*r_gsq))*S_debye(q));

% Cross-sectional form factor
P_CS = @(q) ((2*besselj(1,q*r_cs))/(q*r_cs))^2; 

%------------Calculate-a_long-values-for-correction-function---------------
%NOTE: Translated from c-source code.
%This is hard to decipher what each term means.

p1 = 4.12;
p2 = 4.42;
q0 = 3.1;
c1 = 1.22;
c2 = 0.4288;
c3 = -1.651;
c4 = 1.523;
c5 = 0.1477;
miu = 0.585;
qr_b = q0*r_gsq/b;
qr_b_sq = qr_b^2;
qr_b_4 = qr_b_sq^2;
qr_b_miu = (qr_b)^(-1/miu);
sech2 = (cosh((qr_b - c4)/c5))^(-2);

t1 = (q0^(1+p1+p2))/(b*(p1-p2));
t2 = (c(n_b)/(15*L))*(14*(b^2)*((exp(-qr_b_sq)-1)/(q0*qr_b_sq))+2*q0...
    *r_g*exp(-qr_b_sq)*(11+(7/qr_b_sq)));
t11 = ((c3*qr_b_miu+c2)*qr_b_miu + c1)*qr_b_miu;
t3 = (r_gsq*sech2/(2*c5)*t11);
t4 = r_gsq*(exp(-qr_b_sq)-1+qr_b_sq)*sech2/(c5*qr_b_4);
t5 = (-4*r_gsq*qr_b*(exp(-qr_b_sq)-1))/qr_b_4*(1-w(q0*r_gsq/b));
t10 = 2*(exp(-qr_b_sq)-1+qr_b_sq)/qr_b_4*(1-w(q0*r_gsq/b));
t6 = (4*b/q0)*t10;
t7 = r_gsq*((-3*c3*qr_b_miu-2*c2)*qr_b_miu-c1)*qr_b_miu/(miu*qr_b);
t9 = (c(n_b)/(15*n_b))*(4-exp(-qr_b_sq)*(11+(7/qr_b_sq))+7/qr_b_sq);
t12 = (((b^2)*pi)/(L*q0^2)) +t2+t3-t4+t5-t6+t7*w(q0*r_gsq/b);
t13 = (-b*pi/(L*q0))+t9+t10+t11*w(q0*r_gsq/b);

a1 = (q0^p1)*t13-t1*(q0^(-p2))*(t12+(b*p1)/q0*t13);
a2 = t1*(q0^(-p1))*(t12+b*p1/q0*t13);

a = @(q) a1*((q*b)^(-p1))+ a2*((q*b)^(-p2))+ pi/(q*L);

%-------------------------Piecewise-Functions------------------------------


del = 1.05;
qrange = data(:,1)'; % Define a q range from the minimum q value 
dq = zeros(1,size(qrange,2)-2); % define the derivative values
%for each point, minus the endpoints.

%construct derivative values after value manipulation.
%Centred Difference, Euler's method.

for i = 1:size(qrange,2)-2
    qval1 = qrange(1,i); 
    qval2 = qrange(1,i+2);
    dq(1,i) = (S_exv(qval2)-S_exv(qval1))/(qval2-qval1); 
end

P_CS = @(r,q) (2*besselj(1,q*r)/(q*r))^2; %Redefine

for h = 1:(size(qrange,2)-1)%compute scattering values.
    qvalue = qrange(1,h); %Select q value
    for i = 1:size(rvals,2)
        if L > 4*b %If the Contour length is 4x greater than the Kuhn length
            if qvalue*b <= 3.1
                if dq(1,h) < 0
                   pvals(i,h) = (S_exv(qvals(1,h))+ (c(n_b)/n_b)*(4/15 + 7/(15*u(qvals(1,h)))-(11/15 + 7/(15*u(qvals(1,h))))*exp(-u(qvals(1,h)))));
                else
                   pvals(i,h) = (S_exv_new(qvals(1,h))+ (c(n_b)/n_b)*(4/15 + 7/(15*u(qvals(1,h)))-(11/15 + 7/(15*u(qvals(1,h))))*exp(-u(qvals(1,h)))));
                end
            else
                pvals(i,h) = (a(qvals(1,h)));
            end
        else
            if qvalue*b <= max(1.9/r_gsq,3)
                pvals(i,h) = (S_debye(qvals(1,h))+ (c(n_b)/n_b)*(4/15 + 7/(15*u(qvals(1,h)))-(11/15 + 7/(15*u(qvals(1,h))))*exp(-u(qvals(1,h)))));
            else
                %use a_short from c code

                %This case does not concern us as this is for WLMs with a 
                %Contour length less than 4x the Kuhn length. This occurs at
                %low concentrations and is not relevant at the present time.
            end
        end
    end
end

if B~=0
    if PD ~= 0
        form = pvals;
        for h = 1:(size(qrange,2)-1)
            for i = 1:size(rvals,2)
            form(i,h)=dist_r(1,i)*P_CS(rvals(1,i),qvals(1,h));
            end
        end
        form = mean(form/sum(dist_r/(size(dist_r,2)))); 
    else
        form = pvals;
        for h = 1:(size(qrange,2)-1)
            form(1,h)= P_CS(r_cs,qvals(1,h));
        end
    end
    bc = @(q) B*(sin(q*r_eff)/(q*r_eff))*exp((-(q.^2)*(sigma.^2)));
    sprism = pvals;
    strucfac = qvals;
    for g = 1:size(pvals,2)
        sprism(1,g) = form(1,g)*(pvals(1,g)/(1+bc(qvals(1,g))*pvals(1,g)));
        strucfac(2,g) = 1/(1+bc(qvals(1,g))*pvals(1,g));
    end
    sprism = volume*1e-4*contrast*A*sprism(1,:);
    sprism(2,:) = sprism(1,:);
    sprism(1,:) = pvals(1,:);
end
pvals = volume*1e-4*contrast*A*pvals;

%---------------------------------figures----------------------------------

figure(1)
subplot(3,1,[1 2])
hold on
box on
% title("Scattering function for wormlike micelles. Data file =  "+filename)
set(gca,'Xscale', 'log','fontweight','bold','Linewidth',2.6,'fontsize',14)
set(gca,'Yscale', 'log');
yticks([10^-2 10^0 10^2 10^4]);
xlim([1e-3, 0.5e0])
ylim([1e-3, 1e5])
xlabel("\bf{q (" + Ang + "^{-1})}");
ylabel("\bf{Intensity (cm^{-1})}");

errorbar(data(:,1),data(:,2),data(:,3),'black','marker','s','linewidth',2,'linestyle','none')
plot(qvals(1,1:112),sprism(2,1:112)+BD,'Color',[0 0.5 1],'Linewidth',4)
legend("Data","Model")

residual = data(:,1);
chisq = 0;
for i = 1:size(residual,1)
    exp = data(i,2);
    err = data(i,3);
    fit = sprism(2,i)+BD;
    residual(i,2) = ((exp)-(fit))/(err);
    chisq = chisq + (exp-fit)^2/err^2;
end
r_chisq = chisq/size(sprism,2);

subplot(3,1,3)
hold on 
box on
set(gca,'Xscale', 'log','fontweight','bold','Linewidth',1.6)
xlim([0.1e-2, 0.9e0])
ylim([-3 3])
xlabel("{\it q} (" + Ang + "^{-1})");
ylabel("Normalised Residuals" );
txt = ("{\bf \chi_R ^2 = }"+r_chisq);
text(0.00115,-0.7,txt)

plot(residual(:,1),residual(:,2),'x')
hold off

sprism = sprism';
