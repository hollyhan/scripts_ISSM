%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
steps=[1];

% Logical variable to decide whether to save the file
saveFile = true;

% Define the filename and directory path
fname = 'ln_viscoelastic_Maxwell_compressible_LT120umv0p5lmv2_evaltimes0_5_1000yr_degmax10000';
fpath = '/Users/kyhan/Desktop/LIA-project/Pre-process/Earthmodel/';
fname_ext = strcat(fname, '.mat');
fpath_full = fullfile(fpath, fname_ext);

% define choices for the Earth model solutions
eval_times = [0 : 5 : 1000]; %  [yr] 
rheol_choice = 0 ; % 0 for maxwell, 1 for simple burger, 2 for EBM
thk_litho = 120.0 ; % lithosphere thickness in km
is_compressible = true; % set 'true' if Earth is compressible, 'false' if incompressible.
lowm_vis_coeff = 2.0 ; % coeffitient 'C' for lower mantle viscosity C * 10^21 Pas. e.g. 20 => lower mantle viscosity will be 2X10^22
uppm_vis_coeff = 0.5 ; % coeffitient 'C' for upper mantle viscosity
uppermant_delta = 1.0 ;% for EBM
lowermant_delta = 3.0 ;% for EBM
tauH = 18.6; % in yr, for EBM. Set between month-100 years. 
degmax_sh = 1e4; % spherical harmonics max degree to solve for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath and load necessary files
addpath('/Users/kyhan/Desktop/Data/PREM/PREMearth20240607102846')
load('/Users/kyhan/Desktop/ISSM_REPO/trunk-jpl/test/Data/hypergeom.mat'); 

% register rheological model choice
if rheol_choice == 0
    rheol_name = 'maxwell';
elseif rheol_choice == 1;
    rheol_name = 'simple burger';
elseif rheol_choice == 2
    rheol_name = 'EBM';
end

% Assign the filename to the miscellaneous name
md.miscellaneous.name = fname;

if any(steps==1)
  
    % initialize model
    md=model();
    
    % calculate year to seconds
    yts=365.25*24*3600; % years to seconds. 
    
    % miscellaneous, cluster, etc {{{ 
    md.cluster=generic('name',oshostname(),'np',10);
    md.miscellaneous.name='love';
    md.groundingline.migration='None';
    md.verbose=verbose('all');
    md.verbose=verbose('1111111111');
    % }}} 
    % solid earth structure and materials {{{ 
    md.materials=materials('litho');
    % taken from 3-Publications/23-Adhikari_2021_GRL/2-calculations_v2/1-LoveNumbers/runme_prem.m steps=3; {{{
    radius_id = 1; 
    
    radius_discontinuity = [10 1221.5 3480 3630 5600 5701 5771 5971 6151 6291 6346.6 6356 6368 6371];
    data = csvread('/Users/kyhan/Desktop/Data/PREM/PREMearth20240607102846/PREM_1s.csv'); 
    
    if (radius_id==1)
	    % radius 1 {{{
	    radius = radius_discontinuity; 
	    % }}}
    elseif (radius_id==2)
	    % radius 2 {{{
	    radius = unique(data(:,1)); % km from center of mass 
	    radius(radius<10) = 10; 
	    radius(radius>1221.5 & radius<3480) = 3480; 
	    radius(radius>6368) = 6371; 
	    radius(radius<6368 & radius>6356) = 6368; 
	    radius(radius<6356 & radius>6346.6) = 6346.6; 
	    radius = unique(radius); 
	    radius = unique([radius' radius_discontinuity]); 
	    % }}}
    elseif (radius_id==3)
	    % radius 3 {{{
	    %r0 = linspace(10,6368,80); 
	    r0 = max(10,6372-logspace(0,log10(6371),100)); 
	    r0(r0>=1221.5 & r0<=3480) = 3480;
	    radius = unique([r0 radius_discontinuity]); 
	    % }}}
    end
    
    % extract PREM data (volume-weighted layer properties!) {{{ 
    depth = max(radius)-radius; 
    period=1; % 1s solutions  
    crust =1; % 1 continent, 0 ocean. 
    num_layer = length(radius)-1; 
    
    for jj=1:num_layer;  % from bottom layer to top 
	    rad(1) = radius(jj);	% bottom interface 
	    rad(2) = radius(jj+1); % top interface 
	    if any(rad(2)==radius_discontinuity)
		    rad(2) = rad(2) - 0.01; % to ensure we sample the right properties! 
	    end
	    dep = 6371-rad; 
	    [vpv0,vph0,vsv0,vsh0,rho0] = prem_table_ar(dep,period,crust); 
	    [vpv,vph,vsv,vsh,rho] = prem_volume_average(dep,period,crust); 
	    vpv_avg(jj) = 3*diff(vpv)/diff((rad/6371).^3); 
	    vsv_avg(jj) = 3*diff(vsv)/diff((rad/6371).^3); 
	    rho_avg(jj) = 3*diff(rho)/diff((rad/6371).^3); 
    end
    % }}} 
    
    % plot {{{ 
    vpv_stairs = 0*depth; 
    vpv_stairs(2:end)=vpv_avg; 
    vpv_stairs(1) = vpv_avg(1);
    % 
    rho_stairs = 0*depth; 
    rho_stairs(2:end)=rho_avg; 
    rho_stairs(1) = rho_avg(1);
    % 
    vsv_stairs = 0*depth; 
    vsv_stairs(2:end)=vsv_avg; 
    vsv_stairs(1) = vsv_avg(1);
    
    %figure; 
    %subplot(1,3,1); stairs(rho_stairs,-depth,'-b'); grid on; hold on; 
    %	plot(data(:,3),-data(:,2),'-r'); hold off; 
    %subplot(1,3,2); stairs(vpv_stairs,-depth,'-b'); grid on; hold on; 
    %	plot(data(:,4),-data(:,2),'-r'); hold off; 
    %subplot(1,3,3); stairs(vsv_stairs,-depth,'-b'); grid on; hold on 
    %	plot(data(:,6),-data(:,2),'-r'); hold off; 
    % }}} 
    
    [lambda_avg, mu_avg] = lame(vpv_avg,vsv_avg,rho_avg);
    % }}} 
    
    md.materials.numlayers = num_layer; 
    md.materials.radius = radius'*1e3; % m 
    md.materials.density= rho_avg'*1e3; % kg/m^3
    md.materials.lame_mu= mu_avg'; 
    md.materials.lame_lambda= lambda_avg'; 
    md.materials.issolid = double([(mu_avg'~=0)]); 

    if is_compressible
        str_compressibility = 'compressible';
    else
        % if Earth is incompressible, set lame_lambda parameters to be huge
        md.materials.lame_lambda=md.materials.lame_mu*0+5e17;
        str_compressibility = 'incompressible';
    end
    
    % Maxwell. 
    md.materials.viscosity=ones(md.materials.numlayers,1)*1e21;
    pos = find(radius<3480.0); md.materials.viscosity(pos)=md.materials.viscosity(pos)*0; % core
    pos = find(radius>=3480.0 & radius<(6371.0-670.0)); poslow = pos; % find where lower mantle layer is
    md.materials.viscosity(poslow)=md.materials.viscosity(poslow)*lowm_vis_coeff; % define lower mantle viscosity. 
    pos = find(radius>=(6371.0-670.0) & radius<(6371.0-thk_litho)); posupp = pos; % find where upper mantle layer is
    md.materials.viscosity(posupp)=md.materials.viscosity(posupp)*uppm_vis_coeff; % define upper mantle viscocity.
    pos = find(radius>=(6371.0-thk_litho)); poslitho=pos(1:end-1); % find lithosphere layers.
    md.materials.viscosity(poslitho)=md.materials.viscosity(poslitho)*1e21; % define viscosity lithosphere to be super stiff i.e. elastic
 
    % choose rheology. 
    md.materials.rheologymodel=zeros(md.materials.numlayers,1); % initialize
    md.materials.rheologymodel(:) = rheol_choice; 
    md.materials.rheologymodel(poslitho)=0; % set rheology for lithosphere. Let the lithosphere be maxwell so we dont end up with a small unrelaxed mu for large timescales
    
    % EBM parameters: tauh an delta can have different values for upper and lower mantle. 
    md.materials.ebm_alpha= ones(md.materials.numlayers,1)*0.4; % not so sensitive. ~0.3 probably wont' change '
    md.materials.ebm_delta= ones(md.materials.numlayers,1)*3; % highly sensitive. Between 0.2-5, but open for further exploration. Bigger value => more transient. delta = 0 is Maxwell 
    md.materials.ebm_delta(poslow) = lowermant_delta; 
    md.materials.ebm_delta(posupp) = uppermant_delta; 
    md.materials.ebm_taul= ones(md.materials.numlayers,1)*60*60; %60min*[seconds].Typically between 1 hr - 2 hrs. Dont put it longer than 12 hours. Less effective
    md.materials.ebm_tauh= ones(md.materials.numlayers,1)*tauH*yts; % [seconds] 

    % if large delta, tauh can even go up to 300 for lower mantle., 1-2 years for uppermantle. 
    %md.materials.ebm_tauh(3) = 300.0*yts ;
    %md.materials.ebm_tauh(4) = 1.0*yts ;
    %md.materials.ebm_tauh(5) = 1.0*yts ;

    % Burgers.(just need to exist with the rigth size)
    md.materials.burgers_mu=md.materials.lame_mu/3;
    md.materials.burgers_viscosity=md.materials.viscosity/10;
    % }}} 
    % define love core parameters {{{
    md.love.allow_layer_deletion=1;
    md.love.frequencies=[0];
    md.love.nfreq=length(md.love.frequencies);
    md.love.sh_nmin=1;
    md.love.sh_nmax=degmax_sh;
    md.love.underflow_tol=1e-16;
    md.love.pw_threshold=1e-3;
    md.love.Gravitational_Constant=6.6732e-11;
    md.love.integration_scheme=2; 
    md.love.min_integration_steps=500;
    md.love.max_integration_dr=5e3;
    md.love.allow_layer_deletion=1;
    md.love.forcing_type=11;
    md.love.chandler_wobble=0;
    md.love.complex_computation=0;
    md.love.quad_precision=0;
    
    md.love.istemporal=1;
    md.love.n_temporal_iterations=7; %Prefer 7, and try 6 if you still encounter numerical noise.
    md.love.time = eval_times*yts; %[seconds] 
    
    md.love.hypergeom_table1=h1real;
    md.love.hypergeom_table2=h1real*0;
    md.love.hypergeom_nalpha=size(h1real,1);%101;
    md.love.hypergeom_nz=length(z);
    md.love.hypergeom_z=z;
    
    %fill in all the frequency samples necessary for times requested in md.love.time
    if md.love.istemporal
	    md.love=md.love.build_frequencies_from_time;
    end
    md.love.love_kernels=0;
    % }}} 
    % additional solid earth parameters {{{ 
    md.solidearth.lovenumbers.tk2secular= 0.9668;
    md.solidearth.rotational.equatorialmoi=8.0131e37;
    md.solidearth.rotational.polarmoi=8.0394e37;
    md.solidearth.rotational.angularvelocity=7.292115e-5;
    % }}} 
    md=solve(md,'lv');
    %savemodel(org,md); 
    
    %Loading love numbers
    h = md.results.LoveSolution.LoveHt;
    l = md.results.LoveSolution.LoveLt;
    k = md.results.LoveSolution.LoveKt;
    
    %tidal love numbers (only degree 2 needed for rotational feedback)
    th=md.results.LoveSolution.LoveTidalHt';
    tl=md.results.LoveSolution.LoveTidalLt';
    tk=md.results.LoveSolution.LoveTidalKt';

    %Polar motion transfer function
    pmtf1=md.results.LoveSolution.LovePMTF1t'; %HH this value just for degree 2.  Secular polar motion
    pmtf2=md.results.LoveSolution.LovePMTF2t'; % HH just for chandler wobble
    
    %Time of response (yr) measured since heaviside loading
    time = eval_times; %time in year
    maxwell_time=(md.materials.viscosity./md.materials.lame_mu)/yts ;%HH maxwell time for each layer
    
    % figure; % from Surendra
    % love_d = [2:md.love.sh_nmax]; % 10k. 
    % love_ht =  md.results.LoveSolution.LoveHt(:,3:end); 
    % %love_kt =  1+md.results.LoveSolution.LoveKt(:,3:end); 
    % love_kt_deg =  md.results.LoveSolution.LoveKt(:,3:end).*love_d; 
    % love_lt_deg =  md.results.LoveSolution.LoveLt(:,3:end).*love_d; 
    % subplot(3,1,1),plot(love_d,love_ht); 
    % subplot(3,1,2),plot(love_d,love_lt_deg); 
    % subplot(3,1,3),plot(love_d,love_kt_deg); 
   
    figure;
    %semilogx(time,love_ht(:,3)) % plots degree 2 love number h. HH: h(:,1)=degree 0, h(:,2)=degree 1. 
    semilogx([0:md.love.sh_nmax], h(:,:)); % plot all degrees
    xlabel('n degree','FontSize', 13)
    ylabel('Love Number H','FontSize', 13)
    set(gca, 'FontName', 'Arial', 'FontSize', 12); 
    
    figure;
    semilogx(time, th(3,:), 'DisplayName', 'th','Linewidth', 1)
    hold on
    semilogx(time, tk(3,:), 'DisplayName', 'tk','Linewidth', 1)
    semilogx(time, tl(3,:), 'DisplayName', 'tl','Linewidth', 1)
    legend('show', 'FontSize',13)
    set(gca, 'FontName', 'Arial', 'FontSize', 12); 
    
    figure;
    semilogx([0:md.love.sh_nmax],md.results.LoveSolution(1).LoveHt(1,:),'b','DisplayName','run');
    hold on
    benchmark_prem=lovenumbers('maxdeg',md.love.sh_nmax); % PREM elastic numbers
    semilogx([0:md.love.sh_nmax],benchmark_prem.h,'m', 'DisplayName','benchmark prem');
    legend('show', 'FontSize',13)
    xlabel('n degree','FontSize', 13)
    ylabel('Love Number H','FontSize', 13)
    set(gca, 'FontName', 'Arial', 'FontSize', 12); 
end 

if any(steps==2)
    %%%%%%%%%%% perform the manual picking out of numerically unstable h values in case necessary %%%%%%%%%%%
    h_new=h;
    l_new=l;
    k_new=k;
    
    specific_value =  [120.793; -7.40305; -10.4968; -2.61647; -1.1939; -6.60901; -13.5221; -15.2639]; % change this list depending on the values of unstable solutions
    for i=1:length(specific_value)
        diff = abs(h_new - specific_value(i));    
        [min_difference, linear_index] =min(diff(:));
        [row_index, col_index] = ind2sub(size(h_new), linear_index);
        closest_value = h_new(row_index, col_index);
        h_new(:,col_index)=[]; 
        l_new(:,col_index)=[]; 
        k_new(:,col_index)=[]; 
    end

    figure;
    semilogx(time,h_new(:,:));
    
    % uncomment the below if the new values look satisfying. 
    % h = h_new;
    % l = l_new; 
    % k = k_new;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

ht=h';
lt=l';
kt=k';
tht=ht*0;tht(3,:)=th(3,:);
tkt=kt*0;tkt(3,:)=tk(3,:);
tlt=lt*0;tlt(3,:)=tl(3,:);

%display chosen values 
disp('=============================================')
disp(['rheology model is set to ', rheol_name]);
disp(['Earth is ', str_compressibility]);
fprintf('Lithosphere thickness is set to %i km\n', thk_litho);
fprintf('Upper mantle viscosity is set to %.1e Pas\n', uppm_vis_coeff*1e21);
fprintf('Lower mantle viscosity is set to %.1e Pas\n', lowm_vis_coeff*1e21);
fprintf('Maximum spherical harmonics degree is %i \n', degmax_sh);
disp('=============================================')

if saveFile
    %% save results
    save(fpath_full, 'ht', 'kt', 'lt', 'tht', 'tkt', 'tlt', 'pmtf1', 'pmtf2', 'time');
    disp(['Variables saved to ', fpath_full]);
else
    disp('saveFile is disabled. Not writing into a file.')
end
