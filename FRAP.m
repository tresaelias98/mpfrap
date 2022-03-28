function FRAP(V,g)
% Author: Tresa M. Elias
% this code simulates the molecular dynamics of FRAP by allowing molecules
% to take a random walk on a lattice that is influenced by shear flow. The
% fluorescence recovery curve is then generated and fit to the three
% different FRAP models. 
% Input variables: 
    % V: index of desired velocity (1-27)
    % g: index of desired shear rate (1-48)


tic
% define wz and wr : radial and axial focal volumes, respecitvely 
wz = 5.811*1e-6;
wr = 0.6455*1e-6;

% simulation input variables 
NumMolecules = 20000; % number of desired molecules for simulation
Nreps = 20; % reps
beta = 0.6; % bleach depth parameter
D = 60*1e-12; % input diffusion coefficient
td = (wr^2)/(8*D); % recovery time due to diffusion 

% select particular v and gamma for this simulation using input parameters
vs_input = logspace(-3,2,27); % scaled velocity inputs
gammas_input = logspace(-3,1,48); % scaled shear rate inputs
vs = vs_input(V);
gammas = gammas_input(g);

% convertt scaled velcoity and scaled shear rate to velocity and shear rate
v = vs/(wr/(8*D)); % (vs = vo(wr/8D))
gamma = gammas/((wr^2)/(8*D)); % (ys = y(wr^2/8D))
        
tv = wr/v; % recovery time due to flow
ty = 1/gamma; % recovery time due to shear 
one_over_thalf = 1/td + 1/tv + 1/ty; 
thalf = 1/one_over_thalf; % half recovery time due to diffusion, flow, 
                            % and shear 
tstep = thalf/1000; % define tstep as 1/1000 of the half recovery time 
runtime = thalf*15; % define run time, thalf*15 will ensure the curve is
    % fully recovered witout consuming too much computation time

% determine number of steps each molecules must take
N_step = runtime/tstep; 
% Calculate L (distance each molecules steps)
L = sqrt(6*D*tstep);

% create initial random points 
N_molecules = 0;
count = 1;
    
for k = 1:Nreps
    % while loop runs until there are at least 20000 molecules
    while N_molecules < NumMolecules 
        %10,000 random lattice points 
        factor = 2;
        lattice_x = (-factor*wr:L:factor*wr);
        lattice_y = (-factor*wr:L:factor*wr);
        lattice_z = (-factor*wz:L:factor*wz);

        [X,Y,Z] = meshgrid(lattice_x,lattice_y,lattice_z);
        x = X(:); y = Y(:); z = Z(:);

        % apply probability threshold to each point
        pbl = 1-exp(-beta*exp((-4*(x.^2+y.^2)/(wr^2))-((4*z.^2)/(wz^2))));
        rand_num = rand(length(pbl),1);
        a  = horzcat(x,y,z,pbl,rand_num);
        a = a(a(:,5) < a(:,4),:);

        x1 = a(:,1); y1 = a(:,2); z1 = a(:,3);
        molecules{count} = horzcat(x1,y1,z1);
        N_molecules = N_molecules + length(x1);

        count = count+1;

    end

    % collect all points that passed the probability threshold
    start = vertcat(molecules{1:end});
    % initial start coordinates for all points
    x0 = start(:,1);
    y0 = start(:,2);
    z0 = start(:,3);

    % calculate initial Fbl for this initial distribution
    Fbl = zeros(floor(N_step),1);
    Fbl(1) = sum(exp(((-4*(x0.^2+y0.^2))/(wr^2))+((-4*(z0.^2)/(wz^2)))));
    xnew = zeros([N_molecules 1]);
    ynew = zeros([N_molecules 1]);
    znew = zeros([N_molecules 1]);

    %%%% MAKE MOLECULES TAKE A RANDOM WALK %%%%
    for  m = 2:1:N_step

        % assign each molecule a number between 1 and 6 
        rand_num = randi(6,N_molecules,1);

        % if random integer = 1 : molecule moves in +x
        % [xnew, ynew, znew] = [x0+L, y0, z0]
        znew(rand_num == 1) = z0(rand_num == 1);
        xnew(rand_num == 1) = x0(rand_num == 1) + L + (v+gamma*(znew(rand_num == 1)))*tstep;
        ynew(rand_num == 1) = y0(rand_num == 1);
        % if random integer = 2 : molecule moves in -x
        % [xnew, ynew, znew] = [x0-L, y0, z0]
        znew(rand_num == 2) = z0(rand_num == 2);
        xnew(rand_num == 2) = x0(rand_num == 2) - L + (v+gamma*(znew(rand_num == 2)))*tstep;
        ynew(rand_num == 2) = y0(rand_num == 2);
        % if random integer = 3 : molecule moves in +y
        % [xnew, ynew, znew] = [x0, y0+L, z0]
        znew(rand_num == 3) = z0(rand_num == 3);
        xnew(rand_num == 3) = x0(rand_num == 3) + (v+gamma*(znew(rand_num == 3)))*tstep;
        ynew(rand_num == 3) = y0(rand_num == 3) + L;
        % if random integer = 4 : molecule moves in -y
        % [xnew, ynew, znew] = [x0, y0-L, z0]
        znew(rand_num == 4) = z0(rand_num == 4);
        xnew(rand_num == 4) = x0(rand_num == 4) + (v+gamma*(znew(rand_num == 4)))*tstep;
        ynew(rand_num == 4) = y0(rand_num == 4) - L;
        % if random integer = 5 : molecule moves in +z
        % [xnew, ynew, znew] = [x0, y0, z0+L]
        znew(rand_num == 5) = z0(rand_num == 5) + L;
        xnew(rand_num == 5) = x0(rand_num == 5) + (v+gamma*(znew(rand_num == 5)))*tstep;
        ynew(rand_num == 5) = y0(rand_num == 5);
        % if random integer = 6 : molecule moves in -z
        % [xnew, ynew, znew] = [x0, y0, z0-L]
        znew(rand_num == 6) = z0(rand_num == 6) - L;
        xnew(rand_num == 6) = x0(rand_num == 6) + (v+gamma*(znew(rand_num == 6)))*tstep;
        ynew(rand_num == 6) = y0(rand_num == 6);

        % calculate Fbl
        Fbl_temp = exp(((-4*(xnew.^2+ynew.^2))/(wr^2))+((-4*(znew.^2)/(wz^2))));
        Fbl(m) = sum(Fbl_temp);

        % update new coordinates
        z0 = znew;
        x0 = xnew; 
        y0 = ynew;

    end

    % Calculate F0/Fo
    n = 0:1:10;
    F0Fo = sum((((-beta).^n)./(factorial(n)))./((1+n).^(3/2)));

    % Calculate Fo 
    Fo = Fbl(1)/(1-F0Fo);

    % Calculate FtFo
    FtFo = 1-(Fbl/Fo);
    t = linspace(0,tstep*N_step,length(FtFo))';

    % add noise
    maxCounts = 20000;
    counts = FtFo*maxCounts;
    noise = (poissrnd(counts)./counts);
    FtFo_noise = noise.*FtFo;
    
    % save raw data here

    % determine seed value 
    F = FtFo_noise(1);
    beta_seed = 88.7*F^4 - 246*F^3 + 256*F^2 - 124*F + 25.5; % beta seed value
    half_f = (1+FtFo(1))/2;
    [~,idx] = min(abs(FtFo-half_f));
    tauH = t(idx); % tauH: half recovery time
    tauD = tauH; % tauD: recovery dominated by pure diffusion
    tauV = tauH/sqrt(0.3625); % tauV: recovery dominated by pure flow 
    tauG = tauH/sqrt(0.145); % tauG: recovery dominated by shear flow
    
    % prepare data and options for lsqcurevfit()
    xdata = t;
    ydata = FtFo_noise;
    options = optimset('algorithm', 'levenberg-marquardt',...
                'TolFun', 10^-7,'MaxFunEvals',8*1e3,'MaxIter',1e3,'Display','off');
    
    %%% fit to diffusion only model %%%
    x0_diffusion = [tauD,beta_seed]; % initial guess seed values
    [xfit_diffusion,resnorm_diffusion,~,exitflag,output] = lsqcurvefit(@diffusion_fit,x0_diffusion,xdata,ydata,[],[],options);
    Dfit_diffusion = ((wr*1e6)^2)/(8*xfit_diffusion(1));
    Bfit_diffusion = xfit_diffusion(2);
    
    Ddiffusion(k) = Dfit_diffusion;
    Bdiffusion(k) = Bfit_diffusion;
    Rdiffusion(k) = resnorm_diffusion;
    
    % save diffusion fit data here
    
    %%% fit to diffusion convection model %%%
    x0_convect = [tauD,tauV,beta_seed]; % initial guess seed values
    [xfit_convect,resnorm_convect,~,exitflag,output] = lsqcurvefit(@convective_fit,x0_convect,xdata,ydata,[],[],options);

    Dfit_convect = ((wr*1e6)^2)/(8*xfit_convect(1));
    vfit_convect = (wr*1e6)/xfit_convect(2);
    Bfit_convect = xfit_convect(3);
    
    Dflow(k) = Dfit_convect;
    Vflow(k) = vfit_convect;
    Bflow(k) = Bfit_convect;
    Rflow(k) = resnorm_convect;
    
    % save convective fit data here 

    %%% fit to shear model %%%
    x0_shear = [tauD,tauG,tauV,beta_seed]; % initial guess seed values
    [xfit_shear,resnorm_shear,~,exitflag,output] = lsqcurvefit(@shear_fit,x0_shear,xdata,ydata,[],[],options);
    
    Dfit_shear = ((wr*1e6)^2)/(8*xfit_shear(1));
    vfit_shear = (wr*1e6)/xfit_shear(2);
    gammafit_shear = 1/xfit_shear(3);
    B_shear = xfit_shear(4);
    
    Dshear(k) = Dfit_shear;
    Vshear(k) = vfit_shear;
    Gshear(k) = gammafit_shear;
    Bshear(k) = B_shear;
    Rshear(k) = resnorm_shear;
    
    % save shear data here

end

% finalize results 
Ddiffusion_final = mean(Ddiffusion);
Bdiffusion_final = mean(Bdiffusion);
Rdiffusion_final = mean(Rdiffusion);

Dflow_final = mean(Dflow);
Vflow_final = mean(Vflow);
Bflow_final = mean(Bflow);
Rflow_final = mean(Rflow);

Dshear_final = mean(Dshear);
Vshear_final = mean(Vshear);
Gshear_final = mean(Gshear);
Bshear_final = mean(Bshear);
Rshear_final = mean(Rshear);

% save all final data here 
toc 
