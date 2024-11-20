clear all
close all

%% Universal constants
c=3*10^8;
f_c=15*10^9;
lambda_c=c/f_c;

F= (lambda_c/3)^3;
Gamma=(2*pi*f_c)/50;
%B=7.5*10^9;
B=6*10^9;
f_max= f_c+B/2;
f_min= f_c-B/2;
K=64; %% Number of subcarriers


include_attenuation=1;
atten_coeff=6;

example_no=0;

if example_no==0

    Ny= 8;% 4;
    Nz= 4;% 8;
    theta_0=pi/6;
    %phi_target_arr=[-theta_0+theta_0/Nz:(2*theta_0)/Nz:theta_0-theta_0/Nz];
    delta_dB=3; %0.6;
else
    Ny=  4;
    Nz= 8;
    theta_0=pi/6;
    %phi_target_arr=[-theta_0+theta_0/Nz:(2*theta_0)/Nz:theta_0-theta_0/Nz];
    delta_dB=0.6;
    
end


N=Ny; %% Number of antennas
ft_arr= f_min+ B/K*([0:1:K-1]) + B/(2*K);


%phi_t_arr=[-5.73/180*pi];
%phi_t_arr=[-32.4/180*pi];
phi_t_arr=[-pi/2:pi/100:pi/2];
%phi_t_arr=[-pi/6];
%dy=lambda_c*0.4167;
phi0=pi/6;
%dy=2*c/(2*sin(phi0)*f_max*f_min)*(f_max-f_min);
%dy=1*c/(2*sin(phi0)*f_max*f_min)*(f_max-f_min);
dy=0.4167*lambda_c;
%dy=0.45*lambda_c;
%ng=sin(phi0)*(f_max+f_min)/(f_max-f_min); %% refractive index
ng=2.5;

p_min_arr= f_min*dy/c*(ng+sin(phi_t_arr));
p_max_arr= f_max*dy/c*(ng+sin(phi_t_arr));


ft_star_arr=zeros(length(phi_t_arr),1);
p_star_arr=zeros(length(phi_t_arr),1);

Gdma_ft_phi_t=zeros(length(phi_t_arr), length(ft_arr));

Gdma_starr_ft_star_phi_t=zeros(length(phi_t_arr),1);
Gdma_starr_fc_phi_t=zeros(length(phi_t_arr),1);
Gdma_pin_fc=zeros(length(phi_t_arr),1);
Gdma_pin_all_ON_fc=zeros(length(phi_t_arr),1);

Gdma_starr_ft_star_phi_t_non_ideal = zeros(length(phi_t_arr),1);
Gdma_starr_fc_phi_t_non_ideal = zeros(length(phi_t_arr),1);
Gdma_pin_fc_non_ideal = zeros(length(phi_t_arr),1);


pin_config_arr=zeros(length(phi_t_arr),N);
for phi_i=1:length(phi_t_arr)
    phi_t=phi_t_arr(phi_i);
    
    for f_i=1:length(ft_arr)
        ft=ft_arr(f_i);
        
        X= pi*ft*dy*(ng+sin(phi_t))/c;
        Gdma_ft_phi_t(phi_i, f_i)= 0.25*(N+ abs(sin(N*X)/sin(X)))^2;
        
    end
    if floor(p_max_arr(phi_i))==ceil(p_min_arr(phi_i))
        p_star_arr(phi_i)=floor(p_max_arr(phi_i));
        %p_star_arr(phi_i)=ceil(p_min_arr(phi_i));
        ft_star_arr(phi_i)= p_star_arr(phi_i)*c/(dy*(ng+sin(phi_t)));
    else
        p_feasible_arr=[p_min_arr(phi_i):(p_max_arr(phi_i)-p_min_arr(phi_i))/500 :p_max_arr(phi_i)   ];
        %[~,p_star_i]=max(sin(pi*N*p_feasible_arr)./(sin(pi*p_feasible_arr)));
        [~,p_star_i]=max(abs(sin(pi*N*p_feasible_arr)./(sin(pi*p_feasible_arr))));
        p_star_arr(phi_i)=p_feasible_arr(p_star_i);
        ft_star_arr(phi_i)= p_star_arr(phi_i)*c/(dy*(ng+sin(phi_t)));
    end
    X_star= pi*ft_star_arr(phi_i)*dy*(ng+sin(phi_t))/c;
    Gdma_starr_ft_star_phi_t(phi_i)= 0.25*(N+ abs(sin(N*X_star)/sin(X_star)))^2;
    X_c=pi*f_c*dy*(ng+sin(phi_t))/c;
    Gdma_starr_fc_phi_t(phi_i)=0.25*(N+ abs(sin(N*X_c)/sin(X_c)))^2;
    
    
    
    
    
    
    %% Compute optimal resonant frequencies for a given target angle and target frequency
    f_star_r_l= sqrt(ft_star_arr(phi_i)^2  + Gamma*ft_star_arr(phi_i)/(2*pi)*tan(pi/4*(1-sign(sin(N*X_star)/sin(X_star))  )+  X_star.*([0:1:N-1]-(N-1)/2 )  ));
    alpha_dma_star_f= generate_alpha_weight_vector(ft_star_arr(phi_i), f_star_r_l, F, Gamma, N, 1);
    f_dma_star_f= alpha_dma_star_f.'/(F*(2*pi*f_star_r_l)/Gamma);
    
    f_c_r_l= sqrt(f_c^2  + Gamma*f_c/(2*pi)*tan(pi/4*(1-sign(sin(N*X_c)/sin(X_c))  )+  X_c.*([0:1:N-1]-(N-1)/2 )  ));
    alpha_dma_c_f= generate_alpha_weight_vector(f_c, f_c_r_l, F, Gamma, N, 1);
    f_dma_c_f= alpha_dma_c_f/(F*(2*pi*f_c)/Gamma);
    
    
    
    %% Comparison with PIN diode beamforming
    
    
    Gdma_pin_config_arr=zeros(2^N,1);
    bin_mat= dec2bin([0:1:2^N-1],N-1);
    for n=1:2^N
       Gdma_pin_config_arr(n)= abs(sum(exp(1j*2*X_c*[0:1:N-1]).*(double(bin_mat(n,:)  )- double('0'))  ))^2;
    end
    [max_Gdma_pin_config, n_max]=max(Gdma_pin_config_arr);
    pin_config_arr(phi_i,:) = double(bin_mat(n_max,:)  )- double('0');
    Gdma_pin_fc(phi_i)=max_Gdma_pin_config;
    Gdma_pin_all_ON_fc(phi_i)=abs(sum(exp(1j*2*X_c*[0:1:N-1]).*ones(N,1).'  ))^2;
    
    if include_attenuation==1
        hi_f_non_ideal= generate_intrinsic_phase_vector(ft_star_arr(phi_i), ng, dy, N, c,  atten_coeff, include_attenuation);
        he_f_non_ideal= generate_extrinsic_phase_vector(ft_star_arr(phi_i), dy, N, c, phi_t);
        h_f_non_ideal=he_f_non_ideal.*hi_f_non_ideal;
        alpha_dma_f= generate_alpha_weight_vector(ft_star_arr(phi_i),   f_star_r_l, F, Gamma, N, 1);
        f_dma_f_star= alpha_dma_f/(F*(2*pi*ft_star_arr(phi_i))/Gamma);
        Gdma_starr_ft_star_phi_t_non_ideal(phi_i) = abs(f_dma_f_star.'*h_f_non_ideal)^2;
        
        hi_f_non_ideal_fc= generate_intrinsic_phase_vector(f_c, ng, dy, N, c,  atten_coeff, include_attenuation);
        he_f_non_ideal_fc= generate_extrinsic_phase_vector(f_c, dy, N, c, phi_t);
        h_f_non_ideal_fc=he_f_non_ideal_fc.*hi_f_non_ideal_fc;
        
        Gdma_starr_fc_phi_t_non_ideal(phi_i) = abs( f_dma_c_f.'*h_f_non_ideal_fc)^2;
        Gdma_pin_fc_non_ideal(phi_i) = abs( pin_config_arr(phi_i,:)*h_f_non_ideal_fc)^2;
    end
    
    %% Actual beam pattern
%     
%     phi_arr=[-pi/2:pi/100:pi/2];
%     phi_arr=[-pi/10];
%     Gdma_1_fc=zeros(length(phi_arr),1);
%     Gdma_2_fc=zeros(length(phi_arr),1);
%     Gdma_3_fc=zeros(length(phi_arr),1);
%     Gdma_1_star=zeros(length(phi_arr),1);
%     Gdma_2_star=zeros(length(phi_arr),1);
%     Gdma_3_star=zeros(length(phi_arr),1);
%     f_tilde_dma_c= exp(1j*(-pi/2*sign(sin(N*X_c)/sin(X_c)) +X_c*([0:1:N-1]-(N-1)/2)*2  )).';
%     f_tilde_dma_star= exp(1j*(-pi/2*sign(sin(N*X_star)/sin(X_star)) +X_star*([0:1:N-1]-(N-1)/2)*2 )).';
%     for i=1:length(phi_arr)
%         phi=phi_arr(i);
%         he_c= generate_intrinsic_phase_vector(f_c, ng, dy, N, c);
%         hi_c= generate_extrinsic_phase_vector(f_c, dy, N, c, phi);
%         h_c=he_c.*hi_c;
%         Gdma_1_fc(i)=0.25*abs(-1j*ones(N,1).'*h_c)^2;
%         Gdma_2_fc(i)=0.25*abs(f_tilde_dma_c.'*h_c)^2;
%         Gdma_3_fc(i)= 0.5*real(conj(-1j*ones(N,1).'*h_c)*(f_tilde_dma_c.'*h_c));
%         
%         
%         he_star= generate_intrinsic_phase_vector(ft_star_arr(phi_i), ng, dy, N, c);
%         hi_star= generate_extrinsic_phase_vector(ft_star_arr(phi_i), dy, N, c, phi);
%         h_star=he_star.*hi_star;
%         Gdma_1_star(i)=0.25*abs(-1j*ones(N,1).'*h_star)^2;
%         Gdma_2_star(i)=0.25*abs(f_tilde_dma_star.'*h_star)^2;
%         Gdma_3_star(i)= 0.5*real(conj(-1j*ones(N,1).'*h_star)*(f_tilde_dma_star.'*h_star));
%     end
%     figure
%     plot(phi_arr/pi*180, Gdma_1_fc, 'DisplayName', 'G_{dma,1}(\phi,f_c)');
%     hold on
%     plot(phi_arr/pi*180, Gdma_2_fc,  'DisplayName', 'G_{dma,2}(\phi,f_c)');
%     plot(phi_arr/pi*180, Gdma_3_fc,  'DisplayName', 'G_{dma,3}(\phi,f_c)');
%     plot(phi_arr/pi*180, Gdma_1_fc+Gdma_2_fc+Gdma_3_fc,  'DisplayName', 'G_{dma}(\phi,f_c)');
%     legend
%     ylim([-5,70])
%     grid on
%     xlabel("\phi" +" (in degrees)")
%     ylabel('Beamforming gain (linear scale)')
%     
%     
%     figure
%     plot(phi_arr/pi*180, Gdma_1_star, 'DisplayName', 'G_{dma,1}(\phi,f^{\ast})');
%     hold on
%     plot(phi_arr/pi*180, Gdma_2_star, 'DisplayName', 'G_{dma,2}(\phi,f^{\ast})');
%     plot(phi_arr/pi*180, Gdma_3_star, 'DisplayName', 'G_{dma,3}(\phi,f^{\ast})');
%     plot(phi_arr/pi*180,Gdma_1_star+Gdma_2_star+ Gdma_3_star, 'DisplayName', 'G_{dma}(\phi,f^{\ast})');
%     ylim([-5,70])
%     grid on
%     legend
%     xlabel("\phi" +" (in degrees)")
%     ylabel('Beamforming gain (linear scale)')
%     
    
%     %% Generate frequency response around a bandwidth at the target angle
%     Bt=500*10^6;
%     f_sub= Bt/64;
%     f_arr=[ft_star_arr(phi_i)-Bt/2:  f_sub: ft_star_arr(phi_i)+Bt/2];
%     Gdma1_star_phi_t=zeros(length(f_arr),1);
%     Gdma2_star_phi_t=zeros(length(f_arr),1);
%     Gdma3_star_phi_t=zeros(length(f_arr),1);
%     
%     for fi=1:length(f_arr)
%         he_f= generate_intrinsic_phase_vector(f_arr(fi), ng, dy, N, c);
%         hi_f= generate_extrinsic_phase_vector(f_arr(fi), dy, N, c, phi_t);
%         h_f=he_f.*hi_f;
%         
%         alpha_dma_f= generate_alpha_weight_vector(f_arr(fi), ft_star_arr(phi_i)*ones(N,1), F, Gamma, N, 1);
%         f_dma_f= alpha_dma_f/(F*(2*pi*f_arr(fi))/Gamma);
%         f_tilde_dma_f=2*f_dma_f+1j;
%         
%         Gdma1_star_phi_t(fi)=0.25*abs(-1j*ones(N,1).'*h_f)^2;
%         Gdma2_star_phi_t(fi)=0.25*abs(f_tilde_dma_f.'*h_f)^2;
%         Gdma3_star_phi_t(fi)= 0.5*real(conj(-1j*ones(N,1).'*h_f)*(f_tilde_dma_f.'*h_f));
%         
%     end
%     
%     figure
%     plot(f_arr, Gdma1_star_phi_t, 'DisplayName', 'G_{dma,1}(\phi,f^{\ast})');
%     hold on
%     plot(f_arr, Gdma2_star_phi_t, 'DisplayName', 'G_{dma,2}(\phi,f^{\ast})');
%     plot(f_arr, Gdma3_star_phi_t, 'DisplayName', 'G_{dma,3}(\phi,f^{\ast})');
%     plot(f_arr,Gdma1_star_phi_t+Gdma2_star_phi_t+Gdma3_star_phi_t, 'DisplayName', 'G_{dma}(\phi,f^{\ast})');
%     ylim([-5,70])
%     grid on
%     legend
%     xlabel("Frequency")
%     ylabel('Beamforming gain (linear scale)')
end




figure
plot(phi_t_arr/pi*180, p_min_arr, 'LineWidth', 3, 'DisplayName', "p_{{min}}");
hold on
plot(phi_t_arr/pi*180, p_max_arr, 'LineWidth', 3,  'DisplayName', "p_{{max}}");
plot(phi_t_arr/pi*180, p_star_arr, 'o', 'DisplayName', "p^{\ast}");
grid on
xlabel("\phi" +" (in degrees)")
ylabel("p");
legend

figure
plot(phi_t_arr/pi*180, ft_star_arr/10^9, 'LineWidth', 3);
xlabel("\phi" +" (in degrees)")
ylabel("Target frequency (in GHz)");
grid on



figure
plot(phi_t_arr/pi*180, 10*log10(Gdma_starr_ft_star_phi_t), 'LineWidth', 3, 'DisplayName', 'Optimal target frequency');
hold on
plot(phi_t_arr/pi*180, 10*log10(Gdma_starr_fc_phi_t), 'LineWidth', 3, 'DisplayName', 'Center frequency');
plot(phi_t_arr/pi*180, 10*log10( Gdma_pin_fc), 'LineWidth', 3, 'DisplayName', 'Center frequency (binary)'); 
if include_attenuation==1
    plot(phi_t_arr/pi*180, 10*log10(Gdma_starr_ft_star_phi_t_non_ideal), 'LineWidth', 3, 'DisplayName', 'Optimal target frequency (non ideal DMA)');
    plot(phi_t_arr/pi*180, 10*log10(Gdma_starr_fc_phi_t_non_ideal), 'LineWidth', 3, 'DisplayName', 'Center frequency (non ideal DMA)');
    plot(phi_t_arr/pi*180, 10*log10(Gdma_pin_fc_non_ideal), 'LineWidth', 3, 'DisplayName', 'Center frequency (binary and non ideal DMA)');    
end
xlabel("\phi_t" +" (in degrees)");
ylabel("Beamforming gain (in dB)");
legend
grid on

figure
surf(ft_arr/10^9,  phi_t_arr/pi*180,  10*log10(Gdma_ft_phi_t));
colorbar; caxis([10 max(10*log10(Gdma_ft_phi_t),[], 'all')])
ylabel('Angle (in degrees)')
xlabel('Frequency (in GHz)')
title('SNR response in angle and frequency domain')


%% Planar array simulation


% Ny=4;
% Nz=8;
% theta_0=pi/6;
% phi_target_arr=[-theta_0+theta_0/Nz:(2*theta_0)/Nz:theta_0-theta_0/Nz];
% delta_dB=0.6;
[phi_target_arr, omega_N]= compute_phi_target_arr(ng, Ny, theta_0, delta_dB);
%phi_target_arr=-pi/6;
f_opt_target_arr=c./(dy*(ng+sin(phi_target_arr)));
P=length(phi_target_arr);
Q=Nz/P;
F_res_arr=repmat(f_opt_target_arr, Ny*Q, 1);
f_res_arr= F_res_arr(:);

%phi_arr=[-theta_0:pi/100:theta_0];
phi_0=pi/6;
phi_t_arr_actual=[-phi_0: pi/100:  phi_0];
%phi_arr = [-pi/2:pi/100:pi/2];
fk_arr  = f_min+ B/K*([0:1:K-1]) + B/(2*K);
G_dma_planar_arr=zeros(length(phi_t_arr_actual), length(fk_arr));
%% Generate radiation pattern for planar array
for phi_i=1:length(phi_t_arr_actual)
    phi=phi_t_arr_actual(phi_i);
    for fi=1:length(fk_arr)
        f=fk_arr(fi);
        he_f= generate_intrinsic_phase_vector(f, ng, dy, Ny, c);
        hi_f= generate_extrinsic_phase_vector(f, dy, Ny, c, phi);
        h_f=he_f.*hi_f;
        H_planar_array= repmat(h_f, 1, Nz);
        h_planar_array=H_planar_array(:);
        alpha_dma_f_planar= generate_alpha_weight_vector(f, f_res_arr, F, Gamma, Nz*Ny, 1);
        f_dma_f_planar= alpha_dma_f_planar/(F*(2*pi*f)/Gamma);
        G_dma_planar_arr(phi_i, fi)=abs(f_dma_f_planar.'*h_planar_array)^2;
    end
end




figure
surf( phi_t_arr_actual/pi*180, fk_arr/10^9,   10*log10(G_dma_planar_arr).');
colorbar; caxis([0 max(10*log10(G_dma_planar_arr),[], 'all')])
xlabel('Angle (in degrees)')
ylabel('Frequency (in GHz)')
title('Beamforming gain response in angle and frequency domain')

figure
plot( phi_t_arr_actual/pi*180, max( 10*log10(G_dma_planar_arr), [], 2));

[max_val, max_index] = max( 10*log10(G_dma_planar_arr), [], 2);
figure
plot( phi_t_arr_actual/pi*180, max_index)
xlim([-30, 30])

fk_star=fk_arr(max_index);

figure
plot( phi_t_arr_actual/pi*180, fk_star/10^9)
xlim([-30, 30])
ylim([12, 18])

phi_t_hat=asin(c./(dy.*fk_arr(max_index)) -ng);

figure
plot(phi_t_arr_actual/pi*180, phi_t_hat/pi*180);
xlim([-30, 30])

figure
plot(phi_t_arr_actual/pi*180, abs(phi_t_arr_actual-phi_t_hat)/pi*180);
xlim([-30, 30])
grid on
xlabel('Actual angle (in degrees)')
ylabel('Estimated angle error (in degrees)')

normalized_plot=1;

if normalized_plot==0
    G_phit_fkstar= Nz^2*(sin(pi*Ny*(ng+sin(phi_t_arr_actual))./(ng+sin(phi_t_hat))  )./sin(pi*(ng+sin(phi_t_arr_actual))./(ng+sin(phi_t_hat))  )  ).^2;
    figure
    plot(phi_t_arr_actual/pi*180, 10*log10(G_phit_fkstar),'-o', 'LineWidth', 3, 'DisplayName', "N_y="+num2str(Ny)+", N_z="+num2str(Nz)+' (Target angle estimate)');
    hold on
    plot(phi_t_arr/pi*180, 10*log10(Nz^2*Gdma_starr_ft_star_phi_t),'-x' , 'LineWidth', 3, 'DisplayName', "N_y="+num2str(Ny)+", N_z="+num2str(Nz)+ ' (Exact target angle)');
    xlabel("\phi_t" +" (in degrees)");
    ylabel("Beamforming gain (in dB)");
    xlim([-30, 30])
    legend
    grid on
else
    
    G_phit_fkstar= Nz^2*(sin(pi*Ny*(ng+sin(phi_t_arr_actual))./(ng+sin(phi_t_hat))  )./sin(pi*(ng+sin(phi_t_arr_actual))./(ng+sin(phi_t_hat))  )  ).^2;
    max_Gdma_starr_ft_star_phi_t = max(Gdma_starr_ft_star_phi_t);
    figure
    plot(phi_t_arr_actual/pi*180, 10*log10(G_phit_fkstar)-10*log10(Nz^2*max_Gdma_starr_ft_star_phi_t),'-o', 'LineWidth', 3, 'DisplayName', "N_y="+num2str(Ny)+", N_z="+num2str(Nz)+' (Target angle estimate)');
    hold on
    plot(phi_t_arr/pi*180, 10*log10(Nz^2*Gdma_starr_ft_star_phi_t)-10*log10(Nz^2*max_Gdma_starr_ft_star_phi_t) ,'-x' , 'LineWidth', 3, 'DisplayName', "N_y="+num2str(Ny)+", N_z="+num2str(Nz)+ ' (Exact target angle)');
    xlabel("\phi_t" +" (in degrees)");
    ylabel("Beamforming gain (in dB)");
    xlim([-30, 30])
    legend
    grid on
end 

%% rate calculations with beam training
B_arr=[50:25:500]*10^6; %% data bandwidth
Kd=64; %% number of data subcarriers
Power=250*10^(-3);
kb= 1.380649*10^(-23); %% Boltzmann constant
T=290;  %% Temperature in Kelvin
r=500; %% distance between tx and rx
N_0= kb*T;

Tr_arr = [0.25:0.25:6]*10^9;



rate_dma_starr_ft_star_phi_t=zeros(length(phi_t_arr_actual),length(B_arr), length(Tr_arr));
rate_dma_starr_fc_phi_t=zeros(length(phi_t_arr_actual),length(B_arr), length(Tr_arr));
rate_dma_pin_fc=zeros(length(phi_t_arr_actual),length(B_arr), length(Tr_arr));
rate_dma_pin_all_ON_fc=zeros(length(phi_t_arr_actual),length(B_arr), length(Tr_arr));
rate_dma_beam_trained=zeros(length(phi_t_arr_actual),length(B_arr), length(Tr_arr));
rate_ttd=zeros(length(phi_t_arr_actual),length(B_arr), length(Tr_arr));
rate_ps=zeros(length(phi_t_arr_actual),length(B_arr), length(Tr_arr));

rate_dma_starr_ft_star_phi_t_averaged=zeros(length(B_arr), length(Tr_arr));
rate_dma_starr_fc_phi_t_averaged=zeros(length(B_arr), length(Tr_arr));
rate_dma_pin_fc_averaged=zeros(length(B_arr), length(Tr_arr));
rate_dma_pin_all_ON_fc_averaged=zeros(length(B_arr), length(Tr_arr));
rate_dma_beam_trained_averaged=zeros(length(B_arr), length(Tr_arr));
rate_ttd_averaged =zeros(length(B_arr), length(Tr_arr));
rate_ps_averaged  =zeros(length(B_arr), length(Tr_arr));

ft_star_arr=zeros(length(phi_t_arr_actual), length(Tr_arr));
p_star_arr=zeros(length(phi_t_arr_actual), length(Tr_arr));

for tr_i=1:length(Tr_arr)
    
    Tr= Tr_arr(tr_i);
    
    f_min = f_c-Tr/2;
    f_max = f_c+Tr/2;

    p_min_arr= f_min*dy/c*(ng+sin(phi_t_arr_actual));
    p_max_arr= f_max*dy/c*(ng+sin(phi_t_arr_actual));

    for phi_i=1:length(phi_t_arr_actual)
        phi_t_actual=phi_t_arr_actual(phi_i);
        f_center_beam_trained= fk_star(phi_i);
        f_r_l_m_beam_trained=ones(Ny*Nz, 1)*f_center_beam_trained;


        if floor(p_max_arr(phi_i))==ceil(p_min_arr(phi_i))
            p_star_arr(phi_i, tr_i)=floor(p_max_arr(phi_i));
            %p_star_arr(phi_i)=ceil(p_min_arr(phi_i));
            ft_star_arr(phi_i, tr_i)= p_star_arr(phi_i, tr_i)*c/(dy*(ng+sin(phi_t_actual)));
        else
            p_feasible_arr=[p_min_arr(phi_i):(p_max_arr(phi_i)-p_min_arr(phi_i))/500 :p_max_arr(phi_i)   ];
            %[~,p_star_i]=max(sin(pi*N*p_feasible_arr)./(sin(pi*p_feasible_arr)));
            [~,p_star_i]=max(abs(sin(pi*Ny*p_feasible_arr)./(sin(pi*p_feasible_arr))));
            p_star_arr(phi_i, tr_i)=p_feasible_arr(p_star_i);
            ft_star_arr(phi_i, tr_i)= p_star_arr(phi_i, tr_i)*c/(dy*(ng+sin(phi_t_actual)));
        end
        X_star= pi*ft_star_arr(phi_i, tr_i)*dy*(ng+sin(phi_t_actual))/c;
        f_star_r_l= sqrt(ft_star_arr(phi_i, tr_i)^2  + Gamma*ft_star_arr(phi_i, tr_i)/(2*pi)*tan(pi/4*(1-sign(sin(Ny*X_star)/sin(X_star))  )+  X_star.*([0:1:Ny-1]-(Ny-1)/2 )  ));

        f_center_exact_angle_known= ft_star_arr(phi_i, tr_i);
        %f_r_l_m__exact_angle_known=ones(Ny*Nz, 1)*f_center_exact_angle_known;
        f_r_l_m__exact_angle_known=reshape(repmat(f_star_r_l, Nz,1).', Ny*Nz,1);

        X_c=pi*f_c*dy*(ng+sin(phi_t_actual))/c;
        f_c_r_l= sqrt(f_c^2  + Gamma*f_c/(2*pi)*tan(pi/4*(1-sign(sin(Ny*X_c)/sin(X_c))  )+  X_c.*([0:1:Ny-1]-(Ny-1)/2 )  ));
        f_center_fixed=f_c;
        f_r_l_m_fixed_freq_fc=reshape(repmat(f_c_r_l, Nz,1).', Ny*Nz,1);

         %% Comparison with PIN diode beamforming


        Gdma_pin_config_arr=zeros(2^Ny,1);
        bin_mat= dec2bin([0:1:2^Ny-1],Ny-1);
        for n=1:2^Ny
           Gdma_pin_config_arr(n)= abs(sum(exp(1j*2*X_c*[0:1:Ny-1]).*(double(bin_mat(n,:)  )- double('0'))  ))^2;
        end
        [max_Gdma_pin_config, n_max]=max(Gdma_pin_config_arr);
        pin_config_arr = double(bin_mat(n_max,:)  )- double('0');
        pin_config_arr_l_m=reshape(repmat(pin_config_arr, Nz,1).', Ny*Nz,1);
        f_r_l_m_pin=ones(Ny*Nz, 1)*f_center_fixed;

        for bi=1:length(B_arr)
            BW=B_arr(bi);
            [rate_dma_starr_ft_star_phi_t(phi_i, bi, tr_i), rate_ttd(phi_i, bi, tr_i), rate_ps(phi_i, bi, tr_i)]= compute_achievable_rate_LOS_channel(BW, Kd, Ny, Nz, phi_t_actual,  f_r_l_m__exact_angle_known, f_center_exact_angle_known, F, Gamma, ng, dy, c, Power, N_0,  r,0, atten_coeff, include_attenuation);
            [rate_dma_beam_trained(phi_i, bi, tr_i),~,~]=compute_achievable_rate_LOS_channel(BW, Kd, Ny, Nz, phi_t_actual,  f_r_l_m_beam_trained, f_center_beam_trained, F, Gamma, ng, dy, c, Power, N_0,  r, 0, atten_coeff, include_attenuation);
            [rate_dma_starr_fc_phi_t(phi_i, bi, tr_i), ~,~]=compute_achievable_rate_LOS_channel(BW, Kd, Ny, Nz, phi_t_actual, f_r_l_m_fixed_freq_fc, f_center_fixed, F, Gamma, ng, dy, c, Power, N_0,  r, 0, atten_coeff, include_attenuation);
            [rate_dma_pin_fc(phi_i, bi, tr_i), ~,~]=compute_achievable_rate_LOS_channel(BW, Kd, Ny, Nz, phi_t_actual, f_r_l_m_pin, f_center_fixed, F, Gamma, ng, dy, c, Power, N_0,  r, 1, pin_config_arr_l_m, atten_coeff, include_attenuation);
        end

    end
    rate_dma_beam_trained_averaged(:, tr_i)=mean(rate_dma_beam_trained(:,:,tr_i), 1);
    rate_dma_starr_ft_star_phi_t_averaged(:, tr_i)=mean(rate_dma_starr_ft_star_phi_t(:,:,tr_i)  ,1);
    rate_dma_starr_fc_phi_t_averaged(:, tr_i)=mean(rate_dma_starr_fc_phi_t(:,:,tr_i)  ,1);
    rate_dma_pin_fc_averaged(:, tr_i)=mean(rate_dma_pin_fc(:,:,tr_i),1);
    rate_ttd_averaged(:, tr_i)=mean(rate_ttd(:,:,tr_i), 1);
    rate_ps_averaged(:, tr_i)=mean(rate_ps(:,:,tr_i), 1);

%     figure
%     plot(B_arr/10^6,  rate_dma_starr_ft_star_phi_t_averaged(:, tr_i)/10^6, '-o', 'LineWidth', 3, 'DisplayName', 'Proposed: Optimal target frequency tuning')
%     hold on
%     plot(B_arr/10^6,  rate_dma_beam_trained_averaged(:, tr_i)/10^6, '-v', 'LineWidth', 3, 'DisplayName', 'Proposed: Optimal target frequency tuning (based on angle estimate) ')
%     plot(B_arr/10^6,  rate_dma_starr_fc_phi_t_averaged(:, tr_i)/10^6, '-^', 'LineWidth', 3, 'DisplayName', 'Benchmark: Fixed target frequency tuning (continuous weights)')
%     plot(B_arr/10^6,  rate_dma_pin_fc_averaged(:, tr_i)/10^6, '-sq', 'LineWidth', 3, 'DisplayName', 'Benchmark: Fixed target frequency tuning (Binary weights)')
%     plot(B_arr/10^6,  rate_ttd_averaged(:, tr_i)/10^6, '->', 'LineWidth', 3, 'DisplayName', 'Benchmark: True time delay')
%     plot(B_arr/10^6,  rate_ps_averaged(:, tr_i)/10^6, '-<', 'LineWidth', 3, 'DisplayName', 'Benchmark: Phase shifter')
%     legend
%     grid on
%     xlabel('Bandwidth (MHz)');
%     ylabel('Achievable rate (Mbits/sec)')
%     title("Achievable rate as a function of bandwidth")
end

figure
plot(Tr_arr/10^9,  rate_ttd_averaged(3,:)/10^6, 'r->', 'LineWidth', 3, 'DisplayName', 'Benchmark: True time delay')
hold on
plot(Tr_arr/10^9, rate_dma_starr_ft_star_phi_t_averaged(3,:)/10^6, 'g-o', 'LineWidth', 3, 'DisplayName', 'Proposed: Optimal target frequency tuning')
plot(Tr_arr/10^9,  rate_dma_starr_fc_phi_t_averaged(3,:)/10^6, 'b-sq', 'LineWidth', 3, 'DisplayName', 'Benchmark: Fixed target frequency tuning (continuous weights)')

plot(Tr_arr/10^9,  rate_ttd_averaged(7,:)/10^6, 'r-->', 'LineWidth', 3, 'DisplayName', 'Benchmark: True time delay')
plot(Tr_arr/10^9, rate_dma_starr_ft_star_phi_t_averaged(7,:)/10^6, 'g--o', 'LineWidth', 3, 'DisplayName', 'Proposed: Optimal target frequency tuning')
plot(Tr_arr/10^9,  rate_dma_starr_fc_phi_t_averaged(7,:)/10^6, 'b--sq', 'LineWidth', 3, 'DisplayName', 'Benchmark: Fixed target frequency tuning (continuous weights)')

plot(Tr_arr/10^9,  rate_ttd_averaged(11,:)/10^6, 'r-.>', 'LineWidth', 3, 'DisplayName', 'Benchmark: True time delay')
plot(Tr_arr/10^9, rate_dma_starr_ft_star_phi_t_averaged(11,:)/10^6, 'g-.o', 'LineWidth', 3, 'DisplayName', 'Proposed: Optimal target frequency tuning')
plot(Tr_arr/10^9,  rate_dma_starr_fc_phi_t_averaged(11,:)/10^6, 'b-.sq', 'LineWidth', 3, 'DisplayName', 'Benchmark: Fixed target frequency tuning (continuous weights)')


function  he= generate_intrinsic_phase_vector(f, ng, dy, N, c, atten_coeff, include_attenuation)

he= exp(-1j*ng*2*pi*f/c*dy.*[0:1:N-1]).';


if include_attenuation==1
    he=he.*exp(-atten_coeff*dy.*[0:1:N-1]).';

end

end

function  hi= generate_extrinsic_phase_vector(f, dy, N, c, phi)

hi= exp(-1j*2*pi*f/c*dy.*[0:1:N-1]*sin(phi)).';

end

function alpha_weight_vec= generate_alpha_weight_vector(f, f_res, F, Gamma, N, K)
alpha_weight_vec=zeros(N,K);

for i=1:N
    f_0= f_res(i);
    alpha_weight_vec(i,:)= generate_alpha_weight(f,f_0, Gamma, F);
end

end


function alpha_weight= generate_alpha_weight(f,f_0, Gamma, F)

alpha_weight= 2*pi*f.^2*F./(2*pi*f_0.^2- 2*pi*f.^2+1j*Gamma.*f);

end
function [phi_target_arr, omega_N]=compute_phi_target_arr(ng, Ny, theta_0, delta_dB)
%omega_N=0.44/Ny;
omega_N=compute_cutoff(Ny^2/(10^(delta_dB/10)), Ny);

phi_1= asin((1+omega_N)*(ng+sin(-theta_0))   -ng);

phi_target_arr=phi_1;
phi_index=2;
while phi_target_arr(end)<theta_0
    phi_target_i= asin((1+omega_N)/(1-omega_N)*(ng+sin(phi_target_arr(phi_index-1)))   -ng);
    if phi_target_i<theta_0
        phi_target_arr=[phi_target_arr; phi_target_i];
    else
        break;
    end
    phi_index=phi_index+1;
end
phi_target_arr=phi_target_arr.';
end

function x_cut=  compute_cutoff(t, N)

x_iter= [0:0.001:2/N];

err=abs(abs(sin(N*pi*x_iter)./sin(x_iter*pi)).^2- t);

[~, min_ind]=min(err);
x_cut=x_iter(min_ind);

end


function [rate_DMA, rate_TTD, rate_PS]= compute_achievable_rate_LOS_channel(BW, Kd, Ny, Nz, phi,  f_r_l_m, f_center, F, Gamma, ng, dy, c, P, N_0,  r, varactor_or_pin, pin_config_arr_l_m,  atten_coeff, include_attenuation)

f_min=f_center- BW/2;
f_max=f_center+BW/2;
B_sub=(f_max-f_min)/Kd ;
fd_arr= [f_min:  B_sub: f_max];
SNR_DMA=zeros(length(fd_arr),1);
SNR_TTD=zeros(length(fd_arr),1);

SNR_PS=zeros(length(fd_arr),1);


rate_DMA=0;
rate_TTD=0;
rate_PS=0;

if phi>=0
    optimal_delay= reshape(repmat(dy/c*sin(phi).*[0:1:Ny-1], Nz, 1).', Ny*Nz, 1);
else
    optimal_delay= reshape(repmat(dy/c*sin(phi).*([1:1:Ny]- Ny), Nz, 1).', Ny*Nz, 1);
end

f_PS=  reshape(  repmat(exp(1j*2*pi*f_center*dy/c*sin(phi)*[0:1:Ny-1]), Nz, 1).' , Ny*Nz, 1);

for fdi=1:length(fd_arr)
    f=fd_arr(fdi);
    hi_f= generate_intrinsic_phase_vector(f, ng, dy, Ny, c,  atten_coeff, include_attenuation);
    he_f= generate_extrinsic_phase_vector(f, dy, Ny, c, phi);
    h_f=he_f.*hi_f;
    H_planar_array= repmat(h_f, 1, Nz);
    h_planar_array=H_planar_array(:);
    
    H_ext_planar_array=repmat(he_f, 1, Nz);
    h_ext_planar_array=H_ext_planar_array(:);
    f_TTD= exp(1j.*2*pi*f*optimal_delay);
    
    alpha_dma_f_planar= generate_alpha_weight_vector(f, f_r_l_m, F, Gamma, Nz*Ny, 1);
    f_dma_f_planar= alpha_dma_f_planar/(F*(2*pi*f)/Gamma);
    if varactor_or_pin==1
        f_dma_f_planar=f_dma_f_planar.*pin_config_arr_l_m;
    end
    SNR_DMA(fdi)=(P/(N_0*BW))*abs(f_dma_f_planar.'*h_planar_array)^2*(c/(4*pi*f*r))^2;
    SNR_TTD(fdi)=(P/(N_0*BW))*abs(f_TTD.'*h_ext_planar_array)^2*(c/(4*pi*f*r))^2;
    SNR_PS(fdi)=(P/(N_0*BW))*abs(f_PS.'*h_ext_planar_array)^2*(c/(4*pi*f*r))^2;
    rate_DMA=rate_DMA+  B_sub*log2(1+SNR_DMA(fdi));
    rate_TTD=rate_TTD+B_sub*log2(1+SNR_TTD(fdi));
    rate_PS=rate_PS+B_sub*log2(1+SNR_PS(fdi));
end


end