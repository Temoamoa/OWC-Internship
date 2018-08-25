function [SS_data,Spect,WelevQC] = getRessourceIsurki(Time,DBcell)
    
    Welev = cell2mat(DBcell);
    
    formatOut = 'yyyy-mm-dd HH:MM:SS';

    Tin = datenum(Time.ini,formatOut);
    Tout= datenum(Time.outi,formatOut);
    
    % Define the Tstep
    diff = (Tout-Tin)*24*60*60;
    Ts = diff/length(Welev);
    
    % Collect data  
    Tide = mean(Welev)-7.3; % Set reference tide to the min of measured tides
    
    t_ma = 9;%4*60*2;
    
    tide = zeros(length(Welev)-t_ma,1);
    for i=1:length(Welev)-t_ma
        tide(i)=mean(Welev(i:i+t_ma));
        Welev_tide(i)= Welev(i) - tide(i);
    end

    mean_tide = tide(round(i/2));
    tide(end-t_ma:end) = mean_tide;
    %WelevQC(end-t_ma:end) = Welev(end-t_ma:end) - mean_tide;
    WelevQC = Welev_tide(1:end);%10:end-t_ma-10);
    
    %plot
    figure;plot(tide(1:length(Welev_tide)-t_ma));
    hold on;
    plot(WelevQC);
    plot(Welev(1:length(Welev_tide)-t_ma));
    hold off;
    legend('tide','welevQC','welev')
        
    % Get SS characteristics
    [ f, Sxx ] = specham( Welev_tide', Ts, 16 )	;  %Sxx magnitude, f fréquence
    [Te, Tp]	= specseastate( f, Sxx ); 	
    Hs = H13T13(Welev_tide',Ts);
    
    SS_data = [Hs, Tp, Te, Tide];
    Spect = [ f, Sxx ];
end


%% Ack Joannes Berque 2018
function [f,Sw] = specham(afull,dt,ndof)
% dt		time step in seconds of time series 
% ndof		number of degrees of freedom desired (even number, 2 per chunk of time series) 
nchunk	= 	ndof / 2                        ; % get 2 degrees of freedom per spectral estimate 
Nfull	= 	length( afull ) 				; % time series total record length 

while mod( (Nfull/nchunk/2) ,1)
    fprintf('\n\t Adapt wave elevation TS')
    afull=[afull;0];
    Nfull	= 	length( afull ); 
end
Sxxtot	= 	zeros( Nfull/nchunk/2 + 1, 1 ) 	;  

for idof = 1:nchunk 
	a	= 	afull(  1 + Nfull/nchunk*(idof-1) :  ( Nfull/nchunk * idof ) ) ; 
	n	= 	length(a)							; % number of records, sub-record length is n*dt 
	a	= 	a .* hamming( length(a) ) 			; 
				
	x 					= 	fft( a ) 									; 
	Sxx( 1 ) 			= 	1/n^2	* 	abs( x(1) 		)^2 			; 
	Sxx( 2:(n/2) )		=   2/n^2 	* 	abs( x(2:(n/2)) ).^2 			; 
	Sxx( n/2 + 1 ) 		= 	1/n^2 	* 	abs( x(n/2+1)	)^2 			; 
	
	f					= 		1 / (n*dt) 	* 	( 0 :  (n/2) )			;  

	Sxxtot	= Sxxtot + Sxx' ; 

end % idof = 1:nchunk 

Sxx	= Sxxtot / nchunk 	; 	% variance in spectral bands (m2) 
Sw	= Sxx 	* 	mean( [ afull - mean(afull) ].^2 ) / sum( Sxx ) *(n*dt)	;	% compensate for variance lost tapering  and convert to spectral density (variance per Hz) 
			 
f	= f'	; 	% convert to column vector for convenience 
end

function [Te, Tp] = specseastate( f, S ) 
%  1. spectral significant wave height, m    
%  2. energy mean period, s
%  3. deepwater, linear Airy wave energy flux, kW/m  
%  4. mean absolute wave period, s   
%  5. peak period, s  

% 2016/04/11 - Updating to accept non-uniform spectral interval (to include MEP spectra) 
Df		= 	mean( diff(f) ) 		; 	% mean spectral interval 
rho		= 	1027					; 	% seawater density in kg/m3 
g		= 	9.81					; 	% gravity acceleration m/s2 

% Spectrum and frequency without 0-Hz component 
a_f = find( f ~= 0 );
f1		= 	f( a_f ) 	; 
S1		= 	S( a_f ) 	; 

% Vector of band widths frequencies -- Df(1) is width of band 1 (for f(1)=0 ) 
% 									 -- Df(n)                  n      f(n)   
Nf		= 	length(f) 				; 
Df 		= 	zeros( length(f1), 1 ) 	; 
Df(1) 	= 	f(2)  - f(1) 			; 
Df(Nf) 	= 	f(Nf) - f(Nf-1) 		; 
for jf = 2 : (Nf-1)  
	Df(jf) 	= 		( f(jf+1) 	- f(jf-1) 	) /2  ; 
end 

Df1 = Df( a_f ) ; 

mm1		= 	sum( S1 ./ f1 .* Df1 ) 				; 	% moment -1 of the spectral distribution of the variance in surface elevation 
m0		= 	sum( S .* Df ) 						; 	% 0th moment of the distribution 
m2		= 	sum( S .* f.^2 .*Df ) 				; 	% 2nd moment 

Hm0		= 	4 * sqrt( m0 ) 						; 	% spectral significant wave height 
Te		= 	mm1 / m0							; 	% energy period 
J 		= 		rho * g^2 						...	% deep water energy flux - should be the same as Eflux  
			/ 	(64 * pi) 						... 
			* 	Hm0^2 * Te						... 
			/ 	1000							; 	% convert to kW/m 
Tm02	= 	sqrt( m0 / m2 ) 					; 	% mean wave period in seconds 
[m,i] 	= 	max( S1 )  							; 	% max and index of max in freq spectrum 
Tp		= 	1 / f1(i) 							; 	% peak wave period 

% seastatevec = [ Hm0 			%  1. spectral significant wave height   
% 				Te              %  2. energy mean period 
% 				J               %  3. deepwater, linear Airy wave energy flux 
% 				Tm02            %  4. mean absolute wave period  
% 				Tp	]; 			%  5. peak period 
end

function  tstats = H13T13( elev, Tsamp )
% tstats = [ Hmean, Tzupm, H13, T13, H110, T110 ]

HTz = twvstat( elev, Tsamp ) ; 

Hmean = mean( HTz(:,2) ) ; 
Tzupm = mean( HTz(:,3) ) ; 

% Sorting of wave heights 
[ s1, k ] = sort( HTz(:,2));%, 'r', 'i' ) ; % r each column, i increasing

n13 = round( length(s1) / 3 ) * 2 ; 
H13 = mean( HTz( k( n13:length(s1) ), 2 ) 	)  ; 
T13 = mean( HTz( k( n13:length(s1) ), 3 ) 	)  ; 
	
n110 = round( length(s1) / 10 ) * 9 ; 	
H110 = mean( HTz( k( n110:length(s1) ), 2 )	)  ; 
T110 = mean( HTz( k( n110:length(s1) ), 3 )	)  ; 

% % Sorting of wave periods 
% [ s1, kT ] = gsort( HTz(:,3), 'r', 'i' ) ; 

% tstats = [ Hmean, Tzupm, H13, T13, H110, T110 ] ; 
tstats = H13; 

end

% 	F5.3 	twvstat 	zero upcrossing period and wave heights from time domain 
function HTz = twvstat( h, Ts ) 
% INPUT 
% h		time series of surface elevation (assumed evenly sampled) 
% Ts		sampling interval (e.g. 0.5 sec at bimep) 
% doplot 	optional. Input doplot=1 to plot the time series, wave heights, and zero upcrossings 
% 
% OUTPUT HTz, matrix with the following columns: 
% 		1. 	wave number in the record 
% 		2. 	wave height for wave number iwv in record 
% 		3. 	zero upcrossing period 
% 		4. 	highest surf. elev during this wave 
% 		5. 	lowest  surf. elev during this wave 
% 		6. 	time of beginning of wave 
% 		7. 	time of ending of wave 
% 		8. 	variance in signal during this wave 
% 
% See also waverun (F7.11) 

% 2017/06/01 slight bug fix 
% 2016/04/12 Adding the variance during this wave as output (see checkbimep1/20.4 for use) 
% 2016/03/07 Based on waverun script. Tested and debugged. 
t	= 	( ( 1:length( h ) )'  -1  ) * Ts 		; 	% time of each sample since beginning of record 

% Time of zero-upcrossings 
h0 	= 	[ h	; 	0	] 		; 
h1 	= 	[ 0	; 	h 	] 		; 

iupx 	= 		find( h1 <0 & h0 >=0 ) 	;  
iupx	= 		iupx - 1 				; 
% iupx	= 		iupx( 1: length(iupx)-1 ) 	; 		% 2016/05/09 as per correction in waverunTz2 
													% now no risk of getting spurrious last zero upcross 
tupx	= 		t( iupx ) -	[ 		h( iupx ) 	... 
					./  ( h0( iupx+1 ) - h( iupx ) ) ] * Ts	; 

Nwaves = length( tupx ) - 1 ; 
HTz   = 	[] 	; 	% Matrix with chronological series of each wave hegith and period 
tvec  = 	[] 	; 	% wave times 	for plotting purpose 
hivec = 	[] 	; 	% 				for plotting purpose 
lovec = 	[] 	; 	% 				for plotting purpose 

for iwv = 1 :Nwaves 
	hi 			= 	max( h( iupx(iwv) : iupx(iwv+1) ) ) ; 
	lo 			= 	min( h( iupx(iwv) : iupx(iwv+1) ) ) ; 
	var1		= 	mean( h( iupx(iwv) : iupx(iwv+1) ).^2 ) ; 	% var assuming 0 mean, as in 0 upcrossg 
	HTz	= 	[ 	HTz ; 
				iwv, 							... 	% 1. wave number in the record 
				hi-lo, 							... 	% 2. wave height for wave number iwv in record 
				tupx( iwv+1 ) - tupx( iwv ),  	... 	% 3. zero upcrossing period 
				hi, 							... 	% 4. highest surf. elev during this wave 
				lo								... 	% 5. lowest  surf. elev during this wave 
				tupx( iwv ), 					... 	% 6. time of beginning of wave 
				tupx( iwv+1 ) 					... 	% 7. time of ending of wave 
				var1							] ; 	% 8. variance in signal during this wave 	
end 

end

% F7.0.1 	hamming		returns hamming window for 1D vector
% 2016/01/11. Not sure if the real name is Hanning, Hamming... as these are two different windows
function w = hamming( N ) 

w	= 0.54 - 0.46 * cos( 2 * pi * ( 0 : (N-1) ) / (N-1) )	; 
w = w'; 

end