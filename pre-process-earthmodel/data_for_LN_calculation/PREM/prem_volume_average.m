% PREM model by Victor Tsai (May 3 2012)
%
function [vpv,vph,vsv,vsh,rho] = prem_table_ar(depth,period,crust)
% Input depth in array form (km), single period (s)
vpv = 0*depth;
vph = 0*depth;
vsv = 0*depth;
vsh = 0*depth;
rho = 0*depth;

for i=1:length(depth)
    [vpvt,vpht,vsvt,vsht,rhot]=prem_table(depth(i),period,crust);
    vpv(i) = vpvt;
    vph(i) = vpht;
    vsv(i) = vsvt;
    vsh(i) = vsht;
    rho(i) = rhot;
end

function [vpv1,vph1,vsv1,vsh1,rho] = prem_table(depth,period,crust)
% Input depth in km, period in seconds
% Output vp, vs in km/s; rho in g/cm^3

a=6371;
r = a-depth;
x = r/a;

% Looks up parameters at 1s
if r<0 || r>6371
    error('Depth out of bounds. Exiting.');
elseif r>=0 && r < 1221.5
    rho = 13.0885*x^3/3 - 8.8381*x^5/5;
    vp = 11.2622*x^3/3 - 6.3640*x^5/5;
    vs = 3.6678*x^3/3 - 4.4475*x^5/5;
elseif r>=1221.5 && r<3480
    rho = 12.5815*x^3/3 - 1.2638*x^4/4 - 3.6426*x^5/5 - 5.5281*x^6/6;
    vp = 11.0487*x^3/3 - 4.0362*x^4/4 + 4.8023*x^5/5 - 13.5732*x^6/6;
    vs = 0;
elseif r>=3480 && r<3630
    rho = 7.9565*x^3/3 - 6.4761*x^4/4 + 5.5283*x^5/5 - 3.0807*x^6/6;
    vp = 15.3891*x^3/3 - 5.3181*x^4/4 + 5.5242*x^5/5 - 2.5514*x^6/6;
    vs = 6.9254*x^3/3 + 1.4672*x^4/4 - 2.0834*x^5/5 + 0.9783*x^6/6;
elseif r>=3630 && r<5600
    rho = 7.9565*x^3/3 - 6.4761*x^4/4 + 5.5283*x^5/5 - 3.0807*x^6/6;
    vp = 24.9520*x^3/3 - 40.4673*x^4/4 + 51.4832*x^5/5 - 26.6419*x^6/6;
    vs = 11.1671*x^3/3 - 13.7818*x^4/4 + 17.4575*x^5/5 - 9.2777*x^6/6;
elseif r>=5600 && r<5701
    rho = 7.9565*x^3/3 - 6.4761*x^4/4 + 5.5283*x^5/5 - 3.0807*x^6/6;
    vp = 29.2766*x^3/3 - 23.6027*x^4/4 + 5.5242*x^5/5 - 2.5514*x^6/6;
    vs = 22.3459*x^3/3 - 17.2473*x^4/4 - 2.0834*x^5/5 + 0.9783*x^6/6;
elseif r>=5701 && r<5771
    rho = 5.3197*x^3/3 - 1.4836*x^4/4;
    vp = 19.0957*x^3/3 - 9.8672*x^4/4;
    vs = 9.9839*x^3/3 - 4.9324*x^4/4;
elseif r>=5771 && r<5971
    rho = 11.2494*x^3/3 - 8.0298*x^4/4;
    vp = 39.7027*x^3/3 - 32.6166*x^4/4;
    vs = 22.3512*x^3/3 - 18.5856*x^4/4;
elseif r>=5971 && r<6151
    rho = 7.1089*x^3/3 - 3.8045*x^4/4;
    vp = 20.3926*x^3/3 - 12.2569*x^4/4;
    vs = 8.9496*x^3/3 - 4.4597*x^4/4;
elseif r>=6151 && r<6291
    rho = 2.6910*x^3/3 + 0.6924*x^4/4;
    vpv=0.8317*x^3/3 + 7.2180*x^4/4;
    vph=3.5908*x^3/3 + 4.6172*x^4/4;
    vsv = 5.8582*x^3/3 -1.4678*x^4/4;
    vsh = -1.0839*x^3/3 + 5.7176*x^4/4;
elseif r>=6291 && r<6346.6
    rho = 2.6910*x^3/3 + 0.6924*x^4/4;
    vpv = 0.8317*x^3/3 + 7.2180*x^4/4;
    vph = 3.5908*x^3/3 + 4.6172*x^4/4;
    vsv = 5.8582*x^3/3 - 1.478*x^4/4;
    vsh = -1.0839*x^3/3 + 5.7176*x^4/4;
elseif r>=6346.6 && r<6356
    rho = 2.900*x^3/3;
    vp = 6.800*x^3/3;
    vs = 3.900*x^3/3;
elseif r>=6356 && r<6368
    rho = 2.600*x^3/3;
    vp = 5.800*x^3/3;
    vs = 3.200*x^3/3;
elseif r>=6368 && r<=6371
	 if (crust==1)		% by SA on Nov 12. 2018. 
		rho = 2.600*x^3/3;
		vp = 5.800*x^3/3;
		vs = 3.200*x^3/3;
	 elseif (crust==0)
		rho = 1.020*x^3/3;
		vp = 1.450*x^3/3;
		vs = 0;
	 elseif (crust==0)
		error('Invalid input. what kind of crust? Exiting');
	end
else
    error('Invalid input. Exiting');
end

if r<6151 || r>=6346.6
    vpv = vp;
    vph = vp;
    vsv = vs;
    vsh = vs;
end

% only for 1s. 
vsv1 = vsv; 
vsh1 = vsh; 
vpv1 = vpv; 
vph1 = vph; 


