%==========================================================================
% 2019-2021
% Code written by Ashley P. Williams, SMaCLab, Monash University, Australia

%   ________  ___      ___       __       ______   ___            __       _______   
%  /"       )|"  \    /"  |     /""\     /" _  "\ |"  |          /""\     |   _  "\  
% (:   \___/  \   \  //   |    /    \   (: ( \___)||  |         /    \    (. |_)  :) 
%  \___  \    /\\  \/.    |   /' /\  \   \/ \     |:  |        /' /\  \   |:     \/  
%   __/  \\  |: \.        |  //  __'  \  //  \ _   \  |___    //  __'  \  (|  _  \\  
%  /" \   :) |.  \    /:  | /   /  \\  \(:   _) \ ( \_|:  \  /   /  \\  \ |: |_)  :) 
% (_______/  |___|\__/|___|(___/    \___)\_______) \_______)(___/    \___)(_______/  
                                                                                   

% SMaCLab website can be found here:
% https://sites.google.com/view/smaclab

% Notes:
% Please ensure your file is a 3 column numerical matrix [q, I, I_Error].
% Copy your delimited 3 column data > Workspace drop-down arrow >
% 'Paste' > Set as Numerical Matrix > Rename > Import > Save file in the
% Workspace as the same name in 'Rename' step, in the same directory as this.

% Your final data file, in the form [q, I_model], will be output as
% "PRISMFit" in the 'Workspace' window.

% Your output parameters are found in "PRISM_parameters"
%==========================================================================
clc
clf
clear
Ang = char(197);
%==========================================================================
%                              User Input
%==========================================================================
%==========================================================================
%                    Load File and Make An Initial Guess
%==========================================================================

filename = 'test'; %Insert file name (.mat)
load(filename);
data = test;       %filename (no quotations)

First_point = 2;                % Where to start fitting.
Last_point = size(data, 1);     % Where to stop fitting.

maxEvaluations = 4000;          % Maximum number of fitting evaluations.
toleranceValue = 1e-6;          % Minimum change in function value.

%Input your initial parameter values below

A = 1;         % Scale (Effective Volume Fraction)
r_cs = 20;        % Cross-sectional Radius
axis_ratio = 1;    % Ratio of radii [more info]
BD = 0.001;        % Background.
L = 1000;         % Contour Length
B = 0;           % Peak scale, 0 = no structure factor
b = 100;          % Kuhn Length 
SLD = 1;      % Scattering length density - Cylinder
SLD_solv = 6.3;   % Scattering length density - solvent
sigma = 1;        % Unknown fit parameter for PRISM
r_eff = 60;       % Effective radius

% Dont use PD please. It will break! For a future update.
PD = 0;        % Currently not implemented!

%==========================================================================
%              Manual Fit or Levenberg-Marquardt (L-M) Fit?
%                  1 = L-M fit, 0 = Manual fit
%==========================================================================

LM_fit = 1;

%==========================================================================
%      What would you like to keep fixed? 1 = fixed, 0 = to be fit.
% Parameter order: A, r_cs, axis_ratio, BD, L, B, b, SLD, SLD_solv, sigma, r_eff, PD
%==========================================================================

fixedset = boolean([0,0,1,1,0,0,0,1,1,1,1,1]);

%==========================================================================
%                           End of User Input
%==========================================================================
%==========================================================================
%                  Levenberg-Marquardt (L-M) Fitting
%==========================================================================
q_data = data(First_point:Last_point,1);
I_data = data(First_point:Last_point,2);
Err_data = data(First_point:Last_point,3);

All_parameters = [A, r_cs, axis_ratio, BD, L, B, b, SLD, SLD_solv, sigma, r_eff, PD];
fixedvals = All_parameters(fixedset);
xvary = All_parameters(~fixedset);

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',maxEvaluations,'Tolfun',toleranceValue);
lb = [];
ub = [];
xvary0 = All_parameters(~fixedset);

if LM_fit == 1
    [newParams,resnorm,residual,output] = lsqcurvefit(@(x,data) funwrapper(x, q_data, fixedset,fixedvals),xvary0,q_data,I_data,lb,ub,options);
end

%==========================================================================
%         Plot input and output. Return PRISMfit and PRISM_parameters
%==========================================================================

x(fixedset) = fixedvals;
I_guess = Flexible_cylinder_PRISM(All_parameters, q_data);
PRISM_parameters(:,1) = ["A", "r_cs","axis_ratio", "BD", "L", "B", "b", "SLD", "SLD_solv", "sigma", "r_eff", "PD"];
if LM_fit == 1
    x(~fixedset) = newParams;
    I_model = Flexible_cylinder_PRISM(x, q_data);
    PRISM_parameters(:,2) = x;
else
    I_model = I_guess;
    PRISM_parameters(:,2) = All_parameters;
end
PRISMfit(:,1) = q_data;
PRISMfit(:,2) = I_model;

figure(1)
hold on
box on
set(gca,'Xscale', 'log','fontweight','bold','Linewidth',2.6,'fontsize',14)
set(gca,'Yscale', 'log');
yticks([10^-2 10^0 10^2 10^4]);
xlim([1e-3, 0.7e0])
ylim([1e-3, 1e4])
xlabel("\bf{q (" + Ang + "^{-1})}");
ylabel("\bf{Intensity (cm^{-1})}");

errorbar(q_data(:,1),I_data(:,1),Err_data(:,1),'black','marker','s','linewidth',2,'linestyle','none')
plot(q_data(1:size(q_data)-1),I_guess(1:size(I_guess)-1),'Color',[1 0.5 0],'Linewidth',4)

if LM_fit == 1
    plot(q_data(1:size(q_data)-1),I_model(1:size(I_model)-1),'Color',[0 0.5 1],'Linewidth',4)
    r_chisq = 0;
    chisq = 0;
    for i = 1:size(residual,1)
        chisq = chisq + residual(i,:)^2/Err_data(i,:)^2;
    end
    r_chisq = chisq/(size(q_data,1)-nnz(~fixedset));
end
legend("Data","Initial Guess","L–M fit")

%==========================================================================
%                           Function Wrapper
%==========================================================================

function pred = funwrapper(xvary,q_data,fixedset,fixedvals)
  x = zeros(size(fixedset));
  x(fixedset) = fixedvals;
  x(~fixedset) = xvary;
  pred = Flexible_cylinder_PRISM(x,q_data);
end

%==========================================================================
%             Flexible Cylinder model + PRISM Structure factor
%==========================================================================

function I_model = Flexible_cylinder_PRISM(input, q_data)

A = input(1);
r_cs = input(2);
axis_ratio = input(3);
BD = input(4);
L = input(5);
B = input(6);
b = input(7);
SLD = input(8);
SLD_solv = input(9);
sigma = input(10);
r_eff = input(11);
PD = input(12);


qmin = 0.0001;
qmax = 1;
qvals = q_data(:,1)'; % Range of values
pvals = qvals; 

%==========================================================================
%                       Schulz-Zimm Polydispersity
%==========================================================================
Z = (1-PD^2)/PD^2;
sz = @(x) ((Z+1)^(Z+1))*((x/r_cs)^(Z))*exp(-(Z+1)*(x/r_cs))...
    /(r_cs*gamma(Z+1)); %Schulz distribution function
N_steps = 76;
rvals = ((r_cs-3*PD*r_cs):(6*PD*r_cs)/N_steps:(r_cs+3*PD*r_cs)); %Number of points for PD calculation.
dist_r = [];
for j = 1:size(rvals,2)
    dist_r(1,j) = sz(rvals(1,j)); 
end
dist_r = dist_r/max(dist_r); %Normalise the distribution weights.

%==========================================================================
%                        Initialise other constants
%==========================================================================

contrast = (SLD-SLD_solv).^2; % Del ro^2
volume = pi*(r_cs*axis_ratio)^2*L;

n_b = L./b;

if L > 10*b
    c = @(x) 3.06*(x)^-0.44;
else
    c = @(x) 1;
end

alpha = @(x) sqrt((1+(x/3.12)^2 + (x/8.67)^3)^(0.176/3));

% "Rgsquareshort(L,b),<Rg^2>_0"
r_g0 = (((alpha(n_b)^2)*L*b)/6.0)*(1.0+(b/L)*(-1.5 + (b/L)*(1.5 + (b/L)*0.75*(exp(-2/n_b)-1))));

% "Rgsquare(L,b), <Rg^2>"
r_g = (alpha(n_b)^2)*((b*L)/6);

% Rg = <Rg^2>^0.5
r_gsq = sqrt(r_g);
%r_gsq = sqrt((alpha(n_b)^2)*r_g0);

gauss_N = 76; %According to /lib/Gauss76.c
gauss_W = [.00126779163408536,.00294910295364247,.00462793522803742,.00629918049732845,.00795984747723973,.00960710541471375,.0112381685696677,.0128502838475101,.0144407317482767,.0160068299122486,.0175459372914742,.0190554584671906,.020532847967908,.0219756145344162,.0233813253070112,.0247476099206597,.026072164497986,.0273527555318275,.028587223650054,.029773487255905,.0309095460374916,.0319934843404216,.0330234743977917,.0339977794120564,.0349147564835508,.0357728593807139,.0365706411473296,.0373067565423816,.0379799643084053,.0385891292645067,.0391332242205184,.0396113317090621,.0400226455325968,.040366472122844,.0406422317102947,.0408494593018285,.040987805464794,.0410570369162294,.0410570369162294,.040987805464794,.0408494593018285,.0406422317102947,.040366472122844,.0400226455325968,.0396113317090621,.0391332242205184,.0385891292645067,.0379799643084053,.0373067565423816,.0365706411473296,.0357728593807139,.0349147564835508,.0339977794120564,.0330234743977917,.0319934843404216,.0309095460374916,.029773487255905,.028587223650054,.0273527555318275,.026072164497986,.0247476099206597,.0233813253070112,.0219756145344162,.020532847967908,.0190554584671906,.0175459372914742,.0160068299122486,.0144407317482767,.0128502838475101,.0112381685696677,.00960710541471375,.00795984747723973,.00629918049732845,.00462793522803742,.00294910295364247,.00126779163408536];
gauss_Z = [-.999505948362153,-.997397786355355,-.993608772723527,-.988144453359837,-.981013938975656,-.972229228520377,-.961805126758768,-.949759207710896,-.936111781934811,-.92088586125215,-.904107119545567,-.885803849292083,-.866006913771982,-.844749694983342,-.822068037328975,-.7980001871612,-.77258672828181,-.74587051350361,-.717896592387704,-.688712135277641,-.658366353758143,-.626910417672267,-.594397368836793,-.560882031601237,-.526420920401243,-.491072144462194,-.454895309813726,-.417951418780327,-.380302767117504,-.342012838966962,-.303146199807908,-.263768387584994,-.223945802196474,-.183745593528914,-.143235548227268,-.102483975391227,-.0615595913906112,-.0205314039939986,.0205314039939986,.0615595913906112,.102483975391227,.143235548227268,.183745593528914,.223945802196474,.263768387584994,.303146199807908,.342012838966962,.380302767117504,.417951418780327,.454895309813726,.491072144462194,.526420920401243,.560882031601237,.594397368836793,.626910417672267,.658366353758143,.688712135277641,.717896592387704,.74587051350361,.77258672828181,.7980001871612,.822068037328975,.844749694983342,.866006913771982,.885803849292083,.904107119545567,.92088586125215,.936111781934811,.949759207710896,.961805126758768,.972229228520377,.981013938975656,.988144453359837,.993608772723527,.997397786355355,.999505948362153];
%==========================================================================
%                           Scattering functions
%==========================================================================

%"w_WR(x)
w = @(x) (1+tanh((x-1.523)/0.1477))/2;

u = @(q) r_g*q^2;
u_short = @(q) ((r_g0^0.5)*q)^2;

%Debye scattering
S_debye = @(q) 2*(exp(-u(q))+u(q)-1)/u(q)^2;
S_debye_short = @(q) 2.*(exp(-u_short(q))+u_short(q)-1)/u_short(q)^2;

%Excluded volume scattering
S_exv = @(q) ((1-w(q*r_gsq))*S_debye(q))+(w(q*r_gsq)*(1.22*(q*r_gsq)^(-1/0.585)+0.4288*(q*r_gsq)^(-2/0.585)-1.651*(q*r_gsq)^(-3/0.585)));

%Altered scattering for fixing kink.
S_exv_new = @(q) ((1-w(q*r_gsq))*S_debye(q));

% Cross-sectional form factor
% P_CS = @(q) ((2*besselj(1,q*r_cs))/(q*r_cs))^2; 

%==========================================================================
%              Calculate a_long values for correction function
%==========================================================================
%NOTE: Translated from SasView c-source code. 

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

%==========================================================================
%               Calculate a_short values for correction function
%==========================================================================
%NOTE: Translated from SasView c-source code. 

p1short = 5.36;
p2short = 5.62;
pdiff = p1short-p2short;
r2 = r_g0;
q0short = max(1.9/(r_g0^0.5),3);
exp_qr_b = exp(r2*(q0short/b)^2);
qr2 = q0short*q0short*r2;
b3 = b^3;

q0p = q0short^(-4+p1short);

a1short = ((1.0/(L*r2*r2))*(b/exp_qr_b*q0p)*...
    (8.0*b3*L ...
    - 8.0*b3*exp_qr_b*L ...
    + 2.0*b3*exp_qr_b*L*p2short...
    - 2.0*b*exp_qr_b*L*p2short*qr2 ...
    + 4.0*b*exp_qr_b*L*qr2...
    - 2.0*b3*L*p2short ...
    + 4.0*b*L*qr2 ...
    - pi*exp_qr_b*qr2*q0short*r2...
    + pi*exp_qr_b*p2short*qr2*q0short*r2))/pdiff;

q0p = q0short^(-4+p2short);

a2short = -((1.0/(L*r2*r2))*(b/exp_qr_b*q0p)*...
    (8.0*b3*L ...
    - 8.0*b3*exp_qr_b*L ...
    + 2.0*b3*exp_qr_b*L*p1short...
    - 2.0*b*exp_qr_b*L*p1short*qr2...
    + 4.0*b*exp_qr_b*L*qr2...
    - 2.0*b3*L*p1short...
    + 4.0*b*L*qr2 ...
    - pi*exp_qr_b*qr2*q0short*r2...
    + pi*exp_qr_b*p1short*qr2*q0short*r2))/pdiff;

a_short = @(q) a1short*((q*b)^(-p1short))+ a2short*((q*b)^(-p2short))+ pi/(q*L);

%==========================================================================
%                      Construct piecewise functions
%==========================================================================

del = 1.05;
qrange = q_data(:,1)';          % Define a q range from the minimum q value 
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

for h = 1:(size(qrange,2)-1)     %compute scattering values.
    qvalue = qrange(1,h);        %Select q value
    for i = 1:size(qvals,2)
        if L > 4*b               %If the Contour length is 4x greater than the Kuhn length
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
            if qvalue*b <= q0short
                pvals(i,h) = S_debye_short(qvals(1,h));
            else
                pvals(i,h) =  a_short(qvals(1,h));
            end
        end
    end
end

if PD ~= 0
    %This doesn't work right, maybe modify the rvals first with gaussian,
    %then take schulz polydispersity instead of at the same time.
    form = pvals;
%     valsum = zeros(size(rvals,2),(size(qrange,2)));
%     for h = 1:(size(qrange,2))
%         valsum = zeros(size(qrange,2));
%         for i = 1:size(rvals,2)
%             for k = 1:gauss_N
%                 z = (gauss_Z(k)+1)*(pi/4);
%                 sin_z = sin(z);
%                 cos_z = cos(z);
%                 arg = sqrt(rvals(1,i)^2*sin_z^2 + (rvals(1,i)*axis_ratio)^2*cos_z^2);
%                 valsum(h)=valsum(h) + dist_r(1,i)*gauss_W(k)*P_CS(arg,qvals(1,h))^2;
%             end
%         end
%     end
    %form = mean(0.5.*valsum/sum(dist_r/(size(dist_r,2))));
%     form = mean(form/sum(dist_r/(size(dist_r,2)))); 

    for h = 1:size(qrange,2)
        valsum = 0;
        for k = 1:size(rvals,2)
            for i = 1:gauss_N
                z = (gauss_Z(i)+1)*(pi/4);
                sin_z = sin(z);
                cos_z = cos(z);
                arg = sqrt(rvals(1,k)^2*sin_z^2 + (rvals(1,k)*axis_ratio)^2*cos_z^2);
                valsum = valsum + gauss_W(i)*P_CS(arg,qvals(1,h));
            end
            form(k,h)= 0.5*dist_r(1,k)*valsum;
        end
        form = mean(form/sum(dist_r/(size(dist_r,2)))); 
    end
else
    form = pvals;
    for h = 1:(size(qrange,2))
        valsum = 0;
        for i = 1:gauss_N
            z = (gauss_Z(i)+1)*(pi/4);
            sin_z = sin(z);
            cos_z = cos(z);
            arg = sqrt(r_cs^2*sin_z^2 + (r_cs*axis_ratio)^2*cos_z^2);
            valsum = valsum + gauss_W(i)*P_CS(arg,qvals(1,h));
        end
        form(1,h)= 0.5*valsum;
    end
    
end
if B~=0
    bc = @(q) B*(sin(q*r_eff)/(q*r_eff))*exp((-(q.^2)*(sigma.^2)));
    sprism = pvals;
    strucfac = qvals;
    for g = 1:size(pvals,2)
        sprism(1,g) = form(1,g)*(pvals(1,g)/(1+bc(qvals(1,g))*pvals(1,g)));
        strucfac(2,g) = 1/(1+bc(qvals(1,g))*pvals(1,g));
    end
    sprism = volume*1e-4*contrast*A*sprism;
else
    pvals = volume*1e-4*contrast*A*form(1,:).*pvals;
end

q = qvals(1,:)';
if B~=0
    I_model = sprism(1,:)'+BD;
else
    I_model = pvals(1,:)'+BD;
end
end
%==========================================================================