%% Code to simulate BER for the RAaided LoS MIMO - Ideal and without Quant
clc 
clear all
close all

%% Defining the parameters

lambda = 0.0107; %m
d = 2.84; %m Inter-IRS spacing
fc = 28 * (10^9); %Hz
lambda = 0.0107; %mc
BW = 800 * 10^6;
R = 1.5 * 10^3; %m 
TxDim = 8; %Possible values 8x8 16x16 32x32
RxDim = 8; %Possible values 8x8 16x16 32x32
ds = 0.175; %m
Ntx = 4; % number of tx subarrays
Nrx = 4; % number of rx subarrays

refarr_Msquare = [[],64,256,1024];
num_refarr = length(refarr_Msquare);

z =1; %[[],1:0.1:5]; %m changing different values for depth
num_z = length(z);

M =1024; %The number of RA elements
RA_dim = sqrt(M); %changes based on Number of RA elements for square RA its = sqrt(M)

%% Sending QSPK symbols from the four tx streams and recieving at the Rx
Nsymb = 90000 + 50;
Nsymb_train = 10000;
SISO_SNR_dB_vect = [[],-105:1:80]; %[[],-85:.5:25]
Nsnr = length(SISO_SNR_dB_vect);
SISO_SNR_db_corrected = zeros(Nsnr,1);
Nrun = 10;
BER_MonteCarlo = zeros(Nsnr,Nrun);
BER_MonteCarlo1 = zeros(Nsnr,Nrun);
C_LMMSE = zeros(4,4);
X_est = zeros(4,Nsymb);



%% Setting the locations for the known guys (Tx subarrays)

%Tx-subarray locations using coordinate system (the depth is z)

Tx_array_center = cell(4,1);
Tx_array_center{1} = [ds/2; ds/2; z];
Tx_array_center{2} = [-ds/2; ds/2; z];
Tx_array_center{3} = [-ds/2; -ds/2; z];
Tx_array_center{4} = [ds/2; -ds/2; z];

Tx_array = cell(4,1);
Tx_array_loc = cell(4,1);
base_tx_array = zeros(4,2);

Tx_array{1} = [(ds/2)-(lambda/4),(ds/2)+(lambda/4); (ds/2)+(lambda/4),(ds/2)+(lambda/4); (ds/2)+(lambda/4),(ds/2)-(lambda/4); (ds/2)-(lambda/4),(ds/2)-(lambda/4)]; 
Tx_array{2} = [(-ds/2)-(lambda/4),(ds/2)+(lambda/4); (-ds/2)+(lambda/4),(ds/2)+(lambda/4); (-ds/2)+(lambda/4),(ds/2)-(lambda/4); (-ds/2)-(lambda/4),(ds/2)-(lambda/4)];
Tx_array{3} = [(-ds/2)-(lambda/4),(-ds/2)+(lambda/4); (-ds/2)+(lambda/4),(-ds/2)+(lambda/4); (-ds/2)+(lambda/4),(-ds/2)-(lambda/4); (-ds/2)-(lambda/4),(-ds/2)-(lambda/4)];
Tx_array{4} = [(ds/2)-(lambda/4),(-ds/2)+(lambda/4); (ds/2)+(lambda/4),(-ds/2)+(lambda/4); (ds/2)+(lambda/4),(-ds/2)-(lambda/4); (ds/2)-(lambda/4),(-ds/2)-(lambda/4)];

valc = lambda*0.5; %Inter element spacing
t1 =[];
t2 =[];
t3 = [];
t4 =[];

for iTx_RA = 1:4
    
    tx_loctx =[];
    base_tx_locTX = Tx_array{iTx_RA};
    
    for si = TxDim/2:-1:1
        for sj=TxDim/2:-1:1
            t1 = [t1; base_tx_locTX(1,1)-((sj-1)*valc),base_tx_locTX(1,2)+((si-1)*valc)];
            t2 = [t2; base_tx_locTX(2,1)+((sj-1)*valc),base_tx_locTX(2,2)+((si-1)*valc)];
        end
        tx_loctx = [tx_loctx; t1; flip(t2)];       
        t1=[];
        t2 =[];
    end
    for si = 1:1:TxDim/2
        for sj=TxDim/2:-1:1
            t3 = [t3; base_tx_locTX(4,1)-((sj-1)*valc),base_tx_locTX(4,2)-((si-1)*valc)];
            t4 = [t4; base_tx_locTX(3,1)+((sj-1)*valc),base_tx_locTX(3,2)-((si-1)*valc)];
        end
        tx_loctx = [tx_loctx; t3; flip(t4)];

        t3=[];                       
        t4 =[];
    end
    %%Adding 'depth of z' to the RA loc matrix

    dep_z_tx = zeros(length(tx_loctx),1);
    dep_z_tx(:,1) = z; %Based on the geometry

    tx_loctx = [tx_loctx,dep_z_tx];
    
    Tx_array_loc{iTx_RA} = tx_loctx;

end


%% Finding the steering angle and steering vec at tx and Rx for a particular depth z
R1 = sqrt( ( 1.26^2 ) + (z^2) );
alpha1= asin((1.26)/R1);
theta1 = alpha1; %Elevation angle
phi1 = alpha1; %Azimuth angle

%% Computing the TxDim^2 x 1 Far-field steering vector for the planar subarray
% S_tx = [];
% S_rx = []; 
% for i = 1:TxDim
%     for j = 1:TxDim
%         S_tx =[S_tx; exp(-1j*pi*(((i-1)*sin(theta1)*cos(phi1)) + ((j-1)*sin(theta1)*sin(phi1))))];
%         %S_tx =[S_tx; exp(-1j*pi*(((i-1)*sin(theta1)) + ((j-1)*sin(theta1)*sin(phi1))))];
%     end
% end
% S_tx = S_tx'.';%Taking conjugate since that will be the steering angle
% S_rx = S_tx; %They will be same due to symmetry of the model

%% Computing the TxDim^2 x1 Near-field steering vector for the planar array
% Need to use the centers of the RAs on the tx side basically so lets quickly define that
S_tx = cell(4,1);
S_rx = cell(4,1);

RA_Txside_array_center = cell(4,1);
RA_Txside_array_center{1} = [d/2; d/2; 0];
RA_Txside_array_center{2} = [-d/2; d/2; 0];
RA_Txside_array_center{3} = [-d/2; -d/2; 0];
RA_Txside_array_center{4} = [d/2; -d/2; 0];

for ira = 1:4
    S_txra = [];
    tx_loc = Tx_array_loc{ira};
    RA_loc = RA_Txside_array_center{ira};
    tx_center = Tx_array_center{ira};
    dis_ref = sqrt((tx_center(1) - RA_loc(1))^2 + (tx_center(2) - RA_loc(2))^2 + (tx_center(3) - RA_loc(3))^2);
    for itx =  1:TxDim^2
        dis = sqrt((tx_loc(itx,1) - RA_loc(1))^2 + (tx_loc(itx,2) - RA_loc(2))^2 + (tx_loc(itx,3) - RA_loc(3))^2);
        dis = dis-dis_ref;
        S_txra = [S_txra; exp(-1j*(2*pi/lambda)*dis)];
    end
      % S_tx{ira} = S_txra'.';%Taking conjugate since steering angle
    S_tx{ira} =(1/TxDim)* S_txra'.';%Taking conjugate since steering angle
end
S_rx = S_tx;



for irun = 1:Nrun
%Currently generating a random displacement between -2m to 2m 
RA_disp_x = (-2) + (2-(-2)).*rand(1,1);% The displacement of Rx side RAs and Rx on x axis in m -- horizontal
RA_disp_y =(-2) + (2-(-2)).*rand(1,1); % The displacement of Rx side RAs and Rx on y axis in m -- Vertical


% Setting the location of the recieve subarrays
Rx_array = cell(4,1);
Rx_array_loc = cell(4,1);
base_rx_array = zeros(4,2);

Rx_array{1} = [((ds/2)+RA_disp_x)-(lambda/4),((ds/2)+ RA_disp_y)+(lambda/4); ((ds/2)+RA_disp_x)+(lambda/4),((ds/2)+ RA_disp_y)+(lambda/4); ((ds/2)+RA_disp_x)+(lambda/4),((ds/2)+ RA_disp_y)-(lambda/4); ((ds/2)+RA_disp_x)-(lambda/4),((ds/2)+ RA_disp_y)-(lambda/4)];
Rx_array{2} = [(-(ds/2)+RA_disp_x)-(lambda/4),((ds/2)+ RA_disp_y)+(lambda/4); (-(ds/2)+RA_disp_x)+(lambda/4),((ds/2)+ RA_disp_y)+(lambda/4); (-(ds/2)+RA_disp_x)+(lambda/4),((ds/2)+ RA_disp_y)-(lambda/4); (-(ds/2)+RA_disp_x)-(lambda/4),((ds/2)+ RA_disp_y)-(lambda/4)];
Rx_array{3} = [(-(ds/2)+RA_disp_x)-(lambda/4),(-(ds/2)+ RA_disp_y)+(lambda/4); (-(ds/2)+RA_disp_x)+(lambda/4),(-(ds/2)+ RA_disp_y)+(lambda/4); (-(ds/2)+RA_disp_x)+(lambda/4),(-(ds/2)+ RA_disp_y)-(lambda/4); (-(ds/2)+RA_disp_x)-(lambda/4),(-(ds/2)+ RA_disp_y)-(lambda/4)];
Rx_array{4} = [((ds/2)+RA_disp_x)-(lambda/4),(-(ds/2)+ RA_disp_y)+(lambda/4); ((ds/2)+RA_disp_x)+(lambda/4),(-(ds/2)+ RA_disp_y)+(lambda/4); ((ds/2)+RA_disp_x)+(lambda/4),(-(ds/2)+ RA_disp_y)-(lambda/4); ((ds/2)+RA_disp_x)-(lambda/4),(-(ds/2)+ RA_disp_y)-(lambda/4)];

%Rx center -- we will use this as reference for phase compensation for
%whenr there is vertical or horizontal displacement
Rx_center = [0 + RA_disp_x; 0 + RA_disp_y; R];

% The Rx and Tx subarray "centers" -> useful when considering LoS Channel
% to and from RA and we only need center
Rx_array_center = cell(4,1);
Rx_array_center{1} = [ds/2 + RA_disp_x ; ds/2+ RA_disp_y; R-(z)];
Rx_array_center{2} = [-ds/2+ RA_disp_x; ds/2+ RA_disp_y; R-(z)];
Rx_array_center{3} = [-ds/2+ RA_disp_x; -ds/2+ RA_disp_y; R-(z)];
Rx_array_center{4} = [ds/2+ RA_disp_x; -ds/2+ RA_disp_y; R-(z)];


p1 =[];
p2 =[];
p3 = [];
p4 =[];

for iTx_RA = 1:4
    
    rx_loctx =[];
    base_rx_locTX = Rx_array{iTx_RA};
    
    for si = TxDim/2:-1:1
        for sj=TxDim/2:-1:1
            p1 = [p1; base_rx_locTX(1,1)-((sj-1)*valc),base_rx_locTX(1,2)+((si-1)*valc)];
            p2 = [p2; base_rx_locTX(2,1)+((sj-1)*valc),base_rx_locTX(2,2)+((si-1)*valc)];
        end
        rx_loctx = [rx_loctx; p1; flip(p2)];       
        p1=[];
        p2 =[];
    end
    for si = 1:1:TxDim/2
        for sj=TxDim/2:-1:1
            p3 = [p3; base_rx_locTX(4,1)-((sj-1)*valc),base_rx_locTX(4,2)-((si-1)*valc)];
            p4 = [p4; base_rx_locTX(3,1)+((sj-1)*valc),base_rx_locTX(3,2)-((si-1)*valc)];
        end
        rx_loctx = [rx_loctx; p3; flip(p4)];

        p3=[];                       
        p4 =[];
    end
    %%Adding 'depth of z' to the RA loc matrix

    dep_z_rx = zeros(length(rx_loctx),1);
    dep_z_rx(:,1) = R-z; %Based on the geometry

    rx_loctx = [rx_loctx,dep_z_rx];
    
    Rx_array_loc{iTx_RA} = rx_loctx;

end
t1 =[];
t2 =[];
t3 = [];
t4 =[];
t11 =[];
t22 = [];
t33 =[];
t44 =[];

%% Populating the RAs on the Tx and Rx side

% FYI :Need to change lambda/4 to valc/2 for future like extention to other
% cases of different inter reflect array spacing
%RA Tx side location (The base location for the four RAs in order
%corresponding to the Tx subarrays)
Tx_base_RA_loc = cell(4,1);
Tx_base_RA_loc{1} = [(d/2)-(lambda/4),(d/2)+(lambda/4); (d/2)+(lambda/4),(d/2)+(lambda/4); (d/2)+(lambda/4),(d/2)-(lambda/4); (d/2)-(lambda/4),(d/2)-(lambda/4)];
Tx_base_RA_loc{2} = [-(d/2)-(lambda/4),(d/2)+(lambda/4); -(d/2)+(lambda/4),(d/2)+(lambda/4); -(d/2)+(lambda/4),(d/2)-(lambda/4); -(d/2)-(lambda/4),(d/2)-(lambda/4)];
Tx_base_RA_loc{3} = [-(d/2)-(lambda/4),-(d/2)+(lambda/4); -(d/2)+(lambda/4),-(d/2)+(lambda/4); -(d/2)+(lambda/4),-(d/2)-(lambda/4); -(d/2)-(lambda/4),-(d/2)-(lambda/4)];
Tx_base_RA_loc{4} = [(d/2)-(lambda/4),-(d/2)+(lambda/4); (d/2)+(lambda/4),-(d/2)+(lambda/4); (d/2)+(lambda/4),-(d/2)-(lambda/4); (d/2)-(lambda/4),-(d/2)-(lambda/4)];

%RA Rx side location (I think it will be same as the Tx side and only the z
%value will change, z = R in this case i guess)
Rx_base_RA_loc = cell(4,1);
Rx_base_RA_loc{1} = [((d/2)+RA_disp_x)-(lambda/4),((d/2)+ RA_disp_y)+(lambda/4); ((d/2)+RA_disp_x)+(lambda/4),((d/2)+ RA_disp_y)+(lambda/4); ((d/2)+RA_disp_x)+(lambda/4),((d/2)+ RA_disp_y)-(lambda/4); ((d/2)+RA_disp_x)-(lambda/4),((d/2)+ RA_disp_y)-(lambda/4)];
Rx_base_RA_loc{2} = [(-(d/2)+RA_disp_x)-(lambda/4),((d/2)+ RA_disp_y)+(lambda/4); (-(d/2)+RA_disp_x)+(lambda/4),((d/2)+ RA_disp_y)+(lambda/4); (-(d/2)+RA_disp_x)+(lambda/4),((d/2)+ RA_disp_y)-(lambda/4); (-(d/2)+RA_disp_x)-(lambda/4),((d/2)+ RA_disp_y)-(lambda/4)];
Rx_base_RA_loc{3} = [(-(d/2)+RA_disp_x)-(lambda/4),(-(d/2)+ RA_disp_y)+(lambda/4); (-(d/2)+RA_disp_x)+(lambda/4),(-(d/2)+ RA_disp_y)+(lambda/4); (-(d/2)+RA_disp_x)+(lambda/4),(-(d/2)+ RA_disp_y)-(lambda/4); (-(d/2)+RA_disp_x)-(lambda/4),(-(d/2)+ RA_disp_y)-(lambda/4)];
Rx_base_RA_loc{4} = [((d/2)+RA_disp_x)-(lambda/4),(-(d/2)+ RA_disp_y)+(lambda/4); ((d/2)+RA_disp_x)+(lambda/4),(-(d/2)+ RA_disp_y)+(lambda/4); ((d/2)+RA_disp_x)+(lambda/4),(-(d/2)+ RA_disp_y)-(lambda/4); ((d/2)+RA_disp_x)-(lambda/4),(-(d/2)+ RA_disp_y)-(lambda/4)];

%%Populating the Tx and Rx side RAs with coordinate system 
Tx_RA_loc = cell(4,1);
RX_RA_loc = cell(4,1);
base_RA_locTX = zeros(4,2);
base_RA_locRX = zeros(4,2);

for iTx_RA = 1:4
    
    RA_loctx =[];
    RA_locrx =[];
    base_RA_locTX = Tx_base_RA_loc{iTx_RA};
    base_RA_locRX = Rx_base_RA_loc{iTx_RA};
    
    for si = RA_dim/2:-1:1
        for sj=RA_dim/2:-1:1
            t1 = [t1; base_RA_locTX(1,1)-((sj-1)*valc),base_RA_locTX(1,2)+((si-1)*valc)];
            t2 = [t2; base_RA_locTX(2,1)+((sj-1)*valc),base_RA_locTX(2,2)+((si-1)*valc)];
            %For RX side RAs
            t11 = [t11; base_RA_locRX(1,1)-((sj-1)*valc),base_RA_locRX(1,2)+((si-1)*valc)];
            t22 = [t22; base_RA_locRX(2,1)+((sj-1)*valc),base_RA_locRX(2,2)+((si-1)*valc)];
        end
        RA_loctx = [RA_loctx; t1; flip(t2)];
        RA_locrx = [RA_locrx; t11; flip(t22)];
        t1=[];
        t2 =[];
        t11 = [];
        t22 =[];
    end
    for si = 1:1:RA_dim/2
        for sj=RA_dim/2:-1:1
            t3 = [t3; base_RA_locTX(4,1)-((sj-1)*valc),base_RA_locTX(4,2)-((si-1)*valc)];
            t4 = [t4; base_RA_locTX(3,1)+((sj-1)*valc),base_RA_locTX(3,2)-((si-1)*valc)];
            %For RX side RAs
            t33 = [t33; base_RA_locRX(4,1)-((sj-1)*valc),base_RA_locRX(4,2)-((si-1)*valc)];
            t44 = [t44; base_RA_locRX(3,1)+((sj-1)*valc),base_RA_locRX(3,2)-((si-1)*valc)];
        end
        RA_loctx = [RA_loctx; t3; flip(t4)];
        RA_locrx = [RA_locrx; t33; flip(t44)];
        t3=[];                       
        t4 =[];
        t33 =[];
        t44 = [];
    end
    %%Adding 'depth of z' to the RA loc matrix

    dep_z_tx = zeros(length(RA_loctx),1);
    dep_z_tx(:,1) = 0; %Based on the geometry
    dep_z_rx = zeros(length(RA_locrx),1);
    dep_z_rx(:,1) = R; %Based on the geometry

    RA_loctx = [RA_loctx,dep_z_tx];
    RA_locrx = [RA_locrx,dep_z_rx];
    
    Tx_RA_loc{iTx_RA} = RA_loctx;
    Rx_RA_loc{iTx_RA} = RA_locrx;

end

%% Building the effective 4x4 channel 

% Transmitter side 

H_tx_RA = zeros(4*M, Ntx);
H_tx_RA_phase = zeros(4*M, 4*M);

dis_Tx_RA = cell(4,1);
RA_tx_phase = cell(4,1);
RA_tx_channel = cell(4,1);

%The LoS channel between Tx subarray center to corresponding RA
for idis = 1:4
    dis_mat = zeros(sqrt(M),sqrt(M));
    phase_mat = zeros(sqrt(M),sqrt(M));
    dis_vec = zeros(M,1);
    TxRA_channel = zeros(M,1);
    Tx_loc = Tx_array_center{idis};
    RA_loc = Tx_RA_loc{idis};
    Tx_loc_big = Tx_array_loc{idis};
    
    for di = 1:M
      dis_vec(di) = sqrt( ( ( Tx_loc(1) - RA_loc(di,1))^2 ) + ((Tx_loc(2) - RA_loc(di,2))^2) + ((Tx_loc(3) - RA_loc(di,3))^2));
    end
    temp=1;
    for i = 1:sqrt(M)
        for j= 1:sqrt(M)
            dis_mat(i,j) =  dis_vec(temp,1);
            temp = temp+1;
        end
    end
    
    RA_phase = ((2*pi)/lambda) * dis_mat;
    phase_mat = exp(-1j.*RA_phase);
    dis_Tx_RA{idis} = dis_mat;
    
    temp =1;
    for ic = 1:sqrt(M)
        for jc = 1:sqrt(M)
            TxRA_channel(temp) = ((1/(dis_mat(ic,jc)))*phase_mat(ic,jc))*(2/(4*pi)); %fyi 2 is coming from the root gtgr term 
            temp=temp+1;
        end
    end
    RA_tx_channel{idis} = TxRA_channel;  


% The phase compensation required for the Tx RAs (same needed for the other
% guys on rx side)

%First compensating for phase from Tx subarray to RA


phase_RA_tx1 = zeros(M,1);
dis_txele_ra_mat = zeros(M,TxDim^2);
phase_txele_ra_mat = zeros(M,TxDim^2);



%Using a far-field steering vector 

% for iira = 1:M
%     for iitx = 1:TxDim^2
%         dis_txele_ra_mat(iira,iitx) = sqrt((Tx_loc_big(iitx,1)-RA_loc(iira,1))^2 + (Tx_loc_big(iitx,2)-RA_loc(iira,2))^2 + (Tx_loc_big(iitx,3)-RA_loc(iira,3))^2);
%         phase_txele_ra_mat(iira,iitx) = ((exp(-1j*(2*pi/lambda)*dis_txele_ra_mat(iira,iitx))));
%     end
%     phase_RA_tx1 = phase_txele_ra_mat*((1/TxDim)*S_tx);
% end

%Using near-field steering vector
S_txra = S_tx{idis};
for iira = 1:M
    for iitx = 1:TxDim^2
        dis_txele_ra_mat(iira,iitx) = sqrt((Tx_loc_big(iitx,1)-RA_loc(iira,1))^2 + (Tx_loc_big(iitx,2)-RA_loc(iira,2))^2 + (Tx_loc_big(iitx,3)-RA_loc(iira,3))^2);
        phase_txele_ra_mat(iira,iitx) = ((exp(-1j*(2*pi/lambda)*dis_txele_ra_mat(iira,iitx))));
    end
    phase_RA_tx1 = phase_txele_ra_mat*((1/TxDim)*S_txra);
end





%%Experiment block --> no steering vector, point source
% for iira = 1:M
%         dis_txele_ra = sqrt((Tx_loc(1)-RA_loc(iira,1))^2 + (Tx_loc(2)-RA_loc(iira,2))^2 + (Tx_loc(3)-RA_loc(iira,3))^2);
%         sum_phase =(exp(-1j*(2*pi/lambda)*dis_txele_ra));
%     phase_RA_tx1(iira) = sum_phase; %* sum(S_tx) ;
% end
%%end experiment block

%Calculating the steering angle alpha2 for the compensation towards the rx
%side
tx_rx_center_dis = sqrt((0-Rx_center(1))^2 + (0-Rx_center(2))^2 + (0-Rx_center(3))^2);
alpha2 = acos(R/tx_rx_center_dis);
theta2 = alpha2;
phi2 = alpha2;

S_ra = [];
for i = 1:sqrt(M)
    for j = 1:sqrt(M)
        S_ra =[S_ra; exp(-1j*pi*(((i-1)*sin(theta2)*cos(phi2)) + ((j-1)*sin(theta2)*sin(phi2))))];
    end
end
% S_ra = (1/sqrt(M))*S_ra'.';
S_ra = S_ra'.';

RA_tx_phase{idis} = S_ra.*phase_RA_tx1;


end

%populating the channel matrix from tx to RA


for itx = 1:Ntx
    chan_txRA = RA_tx_channel{itx};
    H_tx_RA((itx-1)*M + 1 : (itx-1)*M + M, itx) = chan_txRA;
end

%Populating the Tx side RA phase matrix
RA_phase_txBig = zeros(4*M,1);
temp=1;
for itx = 1:Ntx
    phase_mat1 = (RA_tx_phase{itx})'.'; %Taking conj as thats what it should compensate for

    for k=1:M
            RA_phase_txBig(temp) = phase_mat1(k);
            temp = temp+1;
    end
end
temp2=1;
for bi = 1:4*M
    for bj = 1:4*M
        if(bi==bj)
           H_tx_RA_phase(bi,bj) =  RA_phase_txBig(temp2);
           temp2 = temp2+1;
        end
    end
end





%Building the Rx-RA channel matrix and the RX side RA phase matrix

% The RX side phase matrix is going to be same as the tx side due to
% symmetry so --> Might change when considering displacement!!!!
H_rx_RA_phase = zeros(4*M, 4*M);
H_rx_RA_phase = H_tx_RA_phase;

RA_rx_phase = cell(4,1);
H_rx_RA = zeros(Ntx,4*M);

for idis = 1:4
    dis_mat = zeros(sqrt(M),sqrt(M));
    dis_vec = zeros(M,1);
    RxRA_channel = zeros(M,1);
    Rx_loc = Rx_array_center{idis};
    RA_loc = Rx_RA_loc{idis};
    phase_mat = zeros(sqrt(M),sqrt(M));
    
    for di = 1:M
      dis_vec(di) = sqrt( ( ( Rx_loc(1) - RA_loc(di,1))^2 ) + ((Rx_loc(2) - RA_loc(di,2))^2) + ((Rx_loc(3) - RA_loc(di,3))^2));       
    end
    temp=1;
    for i = 1:sqrt(M)
        for j= 1:sqrt(M)
            dis_mat(i,j) =  dis_vec(temp,1);
            temp = temp+1;
        end
    end
    
    RA_phase = ((2*pi)/lambda) * dis_mat;
    phase_mat = exp(-1j.*RA_phase);
    
    RA_rx_phase{idis} = phase_mat;
    
    dis_Rx_RA{idis} = dis_mat;

    temp =1;
    for ic = 1:sqrt(M)
        for jc = 1:sqrt(M)
            RxRA_channel(temp) = ((1/(dis_mat(ic,jc)))*phase_mat(ic,jc))*(2/(4*pi));
            temp=temp+1;
        end
    end
    RA_Rx_channel{idis} = RxRA_channel;
    
end

%populating the channel matrix from rx to RA
for itx = 1:Ntx
    chan_rxRA = RA_Rx_channel{itx};
    H_rx_RA(itx,(itx-1)*M + 1 : (itx-1)*M + M) = chan_rxRA;
end

%Building the RA-RA channel (AAAAAHHHHHHH)
H_RA_RA = zeros(4*M,4*M);
disBig = zeros((4*M)^2,1);
phaseBig = zeros((4*M)^2,1);
temp3 = 1;

for itx = 1:Ntx
    tx_RA_pos = Tx_RA_loc{itx}; %isolate RA on tx side
    for i1 =1:M
        %iterating over the M elements from tx side ra
        for irx = 1:4 %isolate RA on Rx side
            rx_RA_pos =  Rx_RA_loc{irx};
            for i2 = 1:M
                %iterating over M elements on Rx side RA
              disBig(temp3) = sqrt( ((tx_RA_pos(i1,1) - rx_RA_pos(i2,1))^2) + ((tx_RA_pos(i1,2) - rx_RA_pos(i2,2))^2) + ((tx_RA_pos(i1,3) - rx_RA_pos(i2,3))^2)); 
              phaseBig(temp3) = exp(-1j*(((2*pi)/lambda) * disBig(temp3)));
              temp3 = temp3+1;
            end
        end
        
    end
    
end

temp4 = 1;
for ira = 1:4*M
    for jra = 1:4*M
        H_RA_RA(ira,jra) = (phaseBig(temp4)/(disBig(temp4)))*(1/(4*pi));
        temp4 = temp4+1;
    end
end

%% The effective 4x4 channel
H_eff = zeros(4,4);
H_eff = H_rx_RA * H_rx_RA_phase * H_RA_RA * H_tx_RA_phase * H_tx_RA;

%% Sending QPSK transmit symbols
Xtx = exp(1i*pi/4 + 1i*pi/2*round(rand(Ntx,Nsymb)*4));
for isnr = 1:Nsnr
        
        X_est = zeros(4,Nsymb);
        SISO_SNR_dB = SISO_SNR_dB_vect(isnr);
        SISO_SNR_lin = 10^(SISO_SNR_dB/10);


        %With Steering vector and subarrays
        %The correct scaling for SISO SNR
        sigmasq = (4*(TxDim^2)^3*(lambda^2))/((4*pi)^2*((R-z)^2)*SISO_SNR_lin); %4 is coming from gain terms
        SISO_SNR_lin_corrected = 1/sigmasq;
        SISO_SNR_db_corrected(isnr) = 10*log10(SISO_SNR_lin_corrected);

        %When you dont use steering vector, use the form below (removing
        %gains from tx and rx subarrays
        % sigmasq = 1/((4*pi)*(R^2)*SISO_SNR_lin);
        % SISO_SNR_lin_corrected = 1/sigmasq;
        % SISO_SNR_db_corrected(isnr) = 10*log10(SISO_SNR_lin_corrected);

        
        AWGN_noise = sqrt(sigmasq/2)*(randn(size(Xtx)) + 1i*randn(size(Xtx)));
        
        Yrx = (H_eff*Xtx) + AWGN_noise;
        %Yrx = Yrx./(3*log2(M^2) + 3*log2((4*M)^2));
        
        for i = 1:4
            C_LMMSE = pinv((H_eff'*H_eff) + (sigmasq*eye(4)) )*H_eff(:,i);
            X_est(i,:) = C_LMMSE'*Yrx;
        end
        
        
        Ph_err_mat = angle(X_est./Xtx);
        BER_LMMSE_MC = (sum(abs(Ph_err_mat)'>pi/4) + sum(abs(Ph_err_mat)'>3*pi/4)*0)/Nsymb/2;
        
        
        %Calculating symbol error rate between the true and estimated
        %values
        
       flag = 0;
       for i = 1:4
           for j = 1:Nsymb
               if sign(real(Xtx(i,j))) == sign(real(X_est(i,j))) && sign(imag(Xtx(i,j))) == sign(imag(X_est(i,j)))
                     flag = flag;
               else
                    flag = flag+1;
               end
           end
       end
       BER_MonteCarlo(isnr,irun) = flag/(4*Nsymb);
       
       
       flag1 = 0;
       for i = 1:4
           for j = 1:Nsymb
               if sign(real(Xtx(i,j))) == sign(real(X_est(i,j))) && sign(imag(Xtx(i,j))) == sign(imag(X_est(i,j)))
                     flag1 = flag1;
               end
               if sign(real(Xtx(i,j))) == sign(real(X_est(i,j))) && sign(imag(Xtx(i,j))) ~= sign(imag(X_est(i,j)))
                    flag1 = flag1+1;
               end
               if sign(real(Xtx(i,j))) ~= sign(real(X_est(i,j))) && sign(imag(Xtx(i,j))) == sign(imag(X_est(i,j)))
                    flag1 = flag1+1;
               end
               if sign(real(Xtx(i,j))) ~= sign(real(X_est(i,j))) && sign(imag(Xtx(i,j))) ~= sign(imag(X_est(i,j)))
                    flag1 = flag1+2;
               end
               
           end
       end
       BER_MonteCarlo1(isnr,irun) = flag1/(4*4*Nsymb);
              
    end
end
ser_vals = zeros(Nsnr,1);
ber_vals = zeros(Nsnr,1);
for i = 1:Nsnr
    ser_vals(i) = sum(BER_MonteCarlo(i,:))/Nrun;
    ber_vals(i) = sum(BER_MonteCarlo1(i,:))/Nrun;
end


% figure
% semilogy(SISO_SNR_dB_vect' + 3*log2(M^2) + 3*log2((4*M)^2), (ser_vals));
% xlabel('Beamformed SNR (dB)')
% ylabel('SER')
% 
% 
% figure
% semilogy(SISO_SNR_dB_vect', (ser_vals));
% xlabel('SNR (dB)')
% ylabel('SER')

% figure
% semilogy(SISO_SNR_dB_vect', (ber_vals));
% xlabel('SNR (dB)')
% ylabel('BER')


figure
semilogy(SISO_SNR_db_corrected, (ber_vals));
xlabel('SNR (dB)')
ylabel('BER')

