
% Modified by Pierre Barrat-Charlaix 10/2017

% Copyright 2014 - by Magnus Ekeberg (magnus.ekeberg@gmail.com)
% All rights reserved
% 
% Permission is granted for anyone to copy, use, or modify this
% software for any uncommercial purposes, provided this copyright 
% notice is retained, and note is made of any changes that have 
% been made. This software is distributed without any warranty, 
% express or implied. In no event shall the author or contributors be 
% liable for any damage arising out of the use of this software.
% 
% The publication of research using this software, modified or not, must include 
% appropriate citations to:
%
% 	M. Ekeberg, C. Lövkvist, Y. Lan, M. Weigt, E. Aurell, Improved contact
% 	prediction in proteins: Using pseudolikelihoods to infer Potts models, Phys. Rev. E 87, 012707 (2013)
%
%	M. Ekeberg, T. Hartonen, E. Aurell, Fast pseudolikelihood
%	maximization for direct-coupling analysis of protein structure
%	from many homologous amino-acid sequences, J. Comput. Phys. 276, 341-356 (2014)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function plmDCA_asymmetric_mask(fastafile_q,outputfile_scores, outputfile_parameters, maskfile,reweighting_threshold,nr_of_cores,lambda,weights_file)
%If should-be numericals are passed as strings, convert them.
    if (isstr(reweighting_threshold))
        reweighting_threshold = str2num(reweighting_threshold);
    end
    if (isstr(nr_of_cores))
        nr_of_cores = str2num(nr_of_cores);
    end

% Minimization options
    options.method='lbfgs'; %Minimization scheme. Default: 'lbfgs', 'cg' for conjugate gradient (use 'cg' if out of RAM).
    options.Display='off';
    options.progTol=1e-3; %Threshold for when to terminate the descent. Default: 1e-9. 
% A note on progTol: In our experiments on PFAM-families, a progTol of 1e-3 gave identical true-positive rates to 1e-9 (default), but with moderately shorter running time. Differences in the scores between progTol 1e-3 and 1e-9 showed up in the 3rd-4th decimal or so (which tends to matter little when ranking them). We here set 1e-7 to be on the safe side, but this can be raised to gain speed. If, however, one wishes to use the scores for some different application, or extract and use the parameters {h,J} directly, we recommend the default progTol 1e-9.

    % addpath(genpath(pwd))

    addpath('/storage1/dgranata/plmDCA_decimation-master/sources/plmDCA_mask/functions')
    addpath('/storage1/dgranata/plmDCA_decimation-master/sources/plmDCA_mask/3rd_party_code/minFunc')
    
% Read inputfile, remove duplicate sequences, and calculate weights and B_eff.
    % Y = dlmread(fastafile_q,' ',1,0)+1;
    Y = dlmread(fastafile_q);
    q = double(max(max(Y)));
    [B_with_id_seq,N]=size(Y);
    Y=unique(Y,'rows');
    [B,N]=size(Y);

% Weights
    weights = ones(B,1);
    if exist('weights_file','var');
        weights = dlmread(weights_file);
        reweighting_threshold = 0;
    end
    if reweighting_threshold>0.0
        fprintf('Starting to calculate weights \n...');
        tic
        %Reweighting in MATLAB:            
        %weights = (1./(1+sum(squareform(pdist(Y,'hamm')<=reweighting_threshold))))';       
         
        %Reweighting in C:
        Y=int32(Y);
        m=calc_inverse_weights(Y-1,reweighting_threshold);
        weights=1./m;

        fprintf('Finished calculating weights \n');
        toc
    end
    B_eff=sum(weights);

    fprintf('### N = %d B_with_id_seq = %d B = %d B_eff = %.2f q = %d\n',N,B_with_id_seq,B,B_eff,q);
   	
% Regularization
    scaled_lambda_h = B_eff*lambda;
    scaled_lambda_J = B_eff/2*lambda;
    fprintf('Effective regularization strength : %f\n',scaled_lambda_h);

% Reading mask    
    try
        masktemp = dlmread(maskfile);
    catch err
        fprintf('mask file was empty, or of incorrect format. No mask will be used\n');
        masktemp = [];
    end
    masksize = size(masktemp,1);
    mask = ones(N,N);
    for s = 1:masksize
        i = masktemp(s,1);
        j = masktemp(s,2);
        mask(i,j) = 0;
        mask(j,i) = 0;
    end
    for i = 1:N
        mask(i,i) = 1;
    end


    Y=int32(Y);q=int32(q);mask = int32(mask);
    w=zeros(q+q^2*(N-1),N); %Matrix in which to store parameter estimates (column r will contain estimates from g_r).
%Run optimizer.
    if nr_of_cores>1
        curr_pool = parpool(nr_of_cores);  
        datetime  
        tic
        parfor r=1:N
            disp(strcat('Minimizing g_r for node r=',int2str(r)))       
            wr=min_g_r(Y,weights,mask(r,:),N,q,scaled_lambda_h,scaled_lambda_J,r,options);
            w(:,r)=wr;
        end
        toc
        delete(curr_pool);
    else
        datetime
        tic
        for r=1:N
            disp(strcat('Minimizing g_r for node r=',int2str(r)))       
            wr=min_g_r(Y,weights,mask(r,:),N,q,scaled_lambda_h,scaled_lambda_J,r,options);
            w(:,r)=wr;
        end
        toc
    end

%Extract the coupling estimates from w.
    % JJ=reshape(w(q+1:end,:),q,q,N-1,N);
    % Jtemp1=zeros(q,q,N*(N-1)/2);
    % Jtemp2=zeros(q,q,N*(N-1)/2);  
    % l=1;
    % for i=1:(N-1)
    %     i
    %      for j=(i+1):N
    %         Jtemp1(:,:,l)=JJ(:,:,j-1,i); %J_ij as estimated from from g_i.
    %         Jtemp2(:,:,l)=JJ(:,:,i,j)'; %J_ij as estimated from from g_j.
    %         l=l+1;
    %     end
    % end

%Shift the coupling estimates into the Ising gauge.
%     J1=zeros(q,q,N*(N-1)/2);
%     J2=zeros(q,q,N*(N-1)/2);
%     for l=1:(N*(N-1)/2)
%         J1(:,:,l)=Jtemp1(:,:,l)-repmat(mean(Jtemp1(:,:,l)),q,1)-repmat(mean(Jtemp1(:,:,l),2),1,q)+mean(mean(Jtemp1(:,:,l)));
%         J2(:,:,l)=Jtemp2(:,:,l)-repmat(mean(Jtemp2(:,:,l)),q,1)-repmat(mean(Jtemp2(:,:,l),2),1,q)+mean(mean(Jtemp2(:,:,l)));
%     end
%     % J1 = Jtemp1;
%     % J2 = Jtemp2;
% %Take J_ij as the average of the estimates from g_i and g_j.
%     J=0.5*(J1+J2);
%     J_NG = 0.5 * (Jtemp1 + Jtemp2);

% %Calculate frob. norms FN_ij.
%     NORMS=zeros(N,N); 
%     l=1;
%     for i=1:(N-1)
%         for j=(i+1):N
%             NORMS(i,j)=norm(J(2:end,2:end,l),'fro');
%             NORMS(j,i)=NORMS(i,j);
%             l=l+1;
%         end
%     end               
       
    
%Calculate scores CN_ij=FN_ij-(FN_i-)(FN_-j)/(FN_--), where '-'
%denotes average.
    % norm_means=mean(NORMS)*N./(sum(mask)-1);
    % % norm_means_all=mean(mean(NORMS))*N/(N-1);
    % norm_means_all = mean(norm_means);
    % CORRNORMS=NORMS-norm_means'*norm_means/norm_means_all;
    % output=[];
    % for i=1:(N-1)
    %     for j=(i+1):N
    %         % output=[output;[i,j,CORRNORMS(i,j)]];
    %         output=[output;[i,j,NORMS(i,j)]];
    %     end
    % end
    % dlmwrite(outputfile_scores,output,'precision',5)
    
    % Parameters in plm gauge
    % Couplings can be set to 0 sum gauge independantly of field values. 
    % Left and right Jij are computed and gauge switched independantly, then averaged in J_tot.
    J_left = zeros(N*q,N*q);
    J_right = zeros(N*q,N*q);
    tic
    for i = 1:N
        for j = (i+1):N
            for a = 1:q
                for b = 1:q
                    J_left((i-1)*q+a, (j-1)*q+b) = w(q + (j-2)*q^2 + (b-1)*q +a,i);
                    J_right((i-1)*q+a, (j-1)*q+b) = w(q + (i-1)*q^2 + (a-1)*q +b,j);
                end
            end
        end
    end
   
    J_left = J_left + J_left';
    J_right = J_right + J_right';
   
    [J_right_0] = switch_gauge(J_right,zeros(1,N*q),'0sum',q);
   
    [J_left_0] = switch_gauge(J_left,zeros(1,N*q),'0sum',q);
   
    J_tot = (J_left_0 + J_right_0)/2;
   

    % h's
    % For each r in 1:N, correponding coupling estimate Jr is computed, then used to switch hr to 0 sum gauge.
    h = zeros(1,N*q);
    Jr = zeros(q,N*q);
    for r = 1:N
        % Constructing couplings estimates for node r
        for i = [1:(r-1),(r+1):N]
            for a = 1:q
                for b = 1:q
                    Jr(a,(i-1)*q+b) = w(q + (i-1-(i>r))*q^2 + (b-1)*q+a,r);
                end
            end
        end
        % Switching gauge of h using those couplings
        h((r-1)*q+(1:q)) = w(1:q,r);
        for i = 1:N
            h((r-1)*q+(1:q)) = h((r-1)*q+(1:q)) + mean(Jr(1:q,(i-1)*q+(1:q)),2)' - mean(mean(Jr(1:q,(i-1)*q+(1:q))));
        end
        h((r-1)*q+(1:q)) = h((r-1)*q+(1:q)) - mean(h((r-1)*q+(1:q)));
    end
% %Calculate frob. norms FN_ij.
     NORMS=zeros(N,N);
     t=zeros(N);
     for i=1:(N-1)
         for j=(i+1):N
             NORMS(i,j)=norm(J_tot((i-1)*q+(2:q),(j-1)*q+(2:q)),'fro');
             NORMS(j,i)=NORMS(i,j);
         end
     end

%Calculate scores CN_ij=FN_ij-(FN_i-)(FN_-j)/(FN_--), where '-'
%denotes average.
     for i=1:(N)
         t(i)=nnz(NORMS(i,:));
     end 
     norm_means=mean(NORMS)*N;
     norm_means_all=mean(mean(NORMS))*N;
     %norm_means_all = mean(norm_means);
     ss=mean(t) %*N
     CORRNORMS=norm_means'*norm_means/norm_means_all;
     output=[];
     for i=1:(N-1)
         for j=(i+1):N
             new_norm=t(i)*t(j)/ss(1)
             output=[output;[i,j,NORMS(i,j)-CORRNORMS(i,j)/new_norm,NORMS(i,j)-CORRNORMS(i,j)/(N-1),NORMS(i,j)]];
             %output=[output;[i,j,NORMS(i,j)]];
         end
     end


    % Writing outputs
    dlmwrite(outputfile_scores,output,'precision',5);
    dlmwrite(outputfile_parameters, [J_tot;h],'delimiter',' ','precision',5);

    toc
end
















function [wr]=min_g_r(Y,weights,mask,N,q,scaled_lambda_h,scaled_lambda_J,r,options)
%Creates function object for (regularized) g_r and minimizes it using minFunc.
    r=int32(r);
    % for i = 1:N
    %     fprintf('mask[%d] = %d\n',i,mask(i));
    % end
    funObj=@(wr)g_r(wr,Y,weights,mask,N,q,scaled_lambda_h,scaled_lambda_J,r);        
    wr0=zeros(q+q^2*(N-1),1);
    wr=minFunc(funObj,wr0,options);    
end

function [fval,grad] = g_r(wr,Y,weights,mask,N,q,lambdah,lambdaJ,r)
%Evaluates (regularized) g_r using the mex-file.
	h_r=reshape(wr(1:q),1,q);
	J_r=reshape(wr(q+1:end),q,q,N-1);

    % if r==1
    %     for a = 1:q
    %         for b = 1:q
    %             fprintf('r = %d -- %f\n',r, J_r(a,b,1));
    %         end
    %     end
    % end

	r=int32(r);
	[fval,grad1,grad2] = g_rC_mask(Y-1,weights,h_r,J_r,[lambdah;lambdaJ],r,mask);
	grad = [grad1(:);grad2(:)];
end

function [N,B,q,Y] = return_alignment(inputfile)
%Reads alignment from inputfile, removes inserts and converts into numbers.
    align_full = fastaread(inputfile);
    B = length(align_full);
    ind = align_full(1).Sequence ~= '.' & align_full(1).Sequence == upper( align_full(1).Sequence );
    N = sum(ind);
    Y = zeros(B,N);

    for i=1:B
        counter = 0;
        for j=1:length(ind)
            if( ind(j) )
                counter = counter + 1;
                Y(i,counter)=letter2number( align_full(i).Sequence(j) );
            end
        end
    end
    q=max(max(Y));
end

function x=letter2number(a)
    switch(a)
        % full AA alphabet
        case '-'
             x=1;
        case 'A'    
            x=2;    
        case 'C'    
            x=3;
        case 'D'
            x=4;
        case 'E'  
            x=5;
        case 'F'
            x=6;
        case 'G'  
            x=7;
        case 'H'
            x=8;
        case 'I'  
            x=9;
        case 'K'
            x=10;
        case 'L'  
            x=11;
        case 'M'
            x=12;
        case 'N'  
            x=13;
        case 'P'
            x=14;
        case 'Q'
            x=15;
        case 'R'
            x=16;
        case 'S'  
            x=17;
        case 'T'
            x=18;
        case 'V'
            x=19;
        case 'W'
            x=20;
        case 'Y'
            x=21;
        otherwise
            x=1;
    end
end












