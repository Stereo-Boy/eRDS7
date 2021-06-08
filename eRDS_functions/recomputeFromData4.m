function psi = recomputeFromData4(psi, plotIt)

%let's recompute everything from all data
psi.new_thresholds = unique([log10([1,10,100,1000,2000,10000]),linspace(log10(0.1),log10(3000),500),log10(psi.history(:,2))']); %range of possible values for thresholds T, in log10(arcsec)              
psi.new_slopes = linspace(0.2,2.4,30);
psi.new_neg_slopes = linspace(0,0.224,20);
trial_nb = size(psi.history,1);
history = nan(trial_nb,16); history(:,1:10)=psi.history;
disparities = psi.new_thresholds;
psi.new_disparities = [1 2]; %not important, will be simplified later
lapse = psi.lapse; g = psi.g; delta = psi.delta; p = psi.p;
tolerance = 10/100; % +-10%

% data structure history
%   1       current psi trial
%   2       current disparity shown in "
%   3       psi.correct or not
%   4-6     current estimates for thres / slope / neg_slope using MAP method
%   7      final 75% threshold estimate using sum parameters (not used)
%   8       trial # (different from psi.trial)
% 9-10      grid range on threshold
% 11-13     pure p(estimates) for thres / slope / neg_slope using MAP method (after marginalizing)
% 14-16     p(estimates +-tolerance) for thres / slope / neg_slope using MAP method (after marginalizing)
[psi.tt, psi.ss, psi.ll, psi.xx] = ndgrid(psi.new_thresholds, psi.new_slopes, psi.new_neg_slopes, psi.new_disparities);
tt=psi.tt(:,:,:,1); ss = psi.ss(:,:,:,1); ll=psi.ll(:,:,:,1);
psi.prior = ones(size(tt))./numel(tt); % reinitialize a flat prior
ttt=tt(:);sss=ss(:);lll=ll(:);
               
                   for i=1:trial_nb % and update with the history of data
                      disp_j = psi.history(i,2); resp_j = psi.history(i,3);
                      psi.likelihoodCR = defineLikelihood_bell(psi.g, ll, ss, psi.delta, psi.p, log10(disp_j), tt, psi.lapse); % psis for success
                      if resp_j == 1
                        % The probability of a response given the disparity shown
                        psi.pCR = sum(psi.likelihoodCR(:).*psi.prior(:)); % success
                        % The posterior, which is the probability of parameters given a response at a potential disparity x
                        psi.prior = psi.likelihoodCR.*psi.prior./psi.pCR; %here?
                      else
                        psi.likelihoodFail = 1 - psi.likelihoodCR; % psis for failure
                        psi.pFail = sum(psi.likelihoodFail(:).*psi.prior(:)); % failure
                        psi.prior = psi.likelihoodFail.*psi.prior./psi.pFail;
                      end
                      [~,idx3] = max(psi.prior(:)); 
                      
                      % max estimate
                       curr_est_max_thr = 10.^ttt(idx3);
                       curr_est_max_pos_slo = sss(idx3);
                       curr_est_max_neg_slo = lll(idx3);
                       history(i,4:6) = [curr_est_max_thr,curr_est_max_pos_slo,curr_est_max_neg_slo];
                       
                       %marginalize distributions for stereoblindness calculation
                        marg_thr=squeeze(sum(sum(psi.prior(:,:,:,1),3),2));  
                        marg_pslo=squeeze(sum(sum(psi.prior(:,:,:,1),3),1)); 
                        marg_nslo=squeeze(sum(sum(psi.prior(:,:,:,1),2),1)); 
                        
                        %calculate and save exact p(MAP)
                        history(i,11:13) = 100.*[marg_thr(disparities==log10(curr_est_max_thr)),marg_pslo(psi.new_slopes==curr_est_max_pos_slo),marg_nslo(psi.new_neg_slopes==curr_est_max_neg_slo)];
                        %calculate and save p(MAP+-5%)
                        selected_thr=marg_thr((disparities>(log10(curr_est_max_thr.*(1-tolerance))))&(disparities<(log10(curr_est_max_thr.*(1+tolerance)))));
                        history(i,14) = 100.*sum(selected_thr(:)); % probability that this threshold is correct (cumulated prob sum of t+-tolerance)
                        selected_slo=marg_pslo((psi.new_slopes>(curr_est_max_pos_slo.*(1-tolerance)))&(psi.new_slopes<(curr_est_max_pos_slo.*(1+tolerance))));
                        history(i,15) = 100.*sum(selected_slo(:));
                        selected_nslo=marg_nslo((psi.new_neg_slopes>(curr_est_max_neg_slo.*(1-tolerance)))&(psi.new_neg_slopes<(curr_est_max_neg_slo.*(1+tolerance))));
                        history(i,16) = 100.*sum(selected_nslo(:));
                       if plotIt==1 && i==trial_nb
                            psi.stereoblind_prob = 100*sum(marg_thr((10.^psi.new_thresholds)>=psi.maxAllowerThreshold));   
                            
                            %calculate CI for threshold from marginalized posterior distribution
                             pp=0; k=1; idx = find(disparities==log10(curr_est_max_thr));
                             while pp<95
                                if (((idx-k)>0) && ((idx+k)<=numel(disparities))); select_p=marg_thr((idx-k):(idx+k)); else; disp('Asymetric CI - discard them'); k=k-1; break; end
                                pp = 100.*sum(select_p);
                                k=k+1;
                             end
                             thr_CI = 10.^[disparities(idx-k),disparities(idx+k)];
                             disp(['Threshold confidence interval (95%) = [',sprintf('%2.f"',thr_CI(1)),'-',sprintf('%2.f"',thr_CI(2)),']'])
                                subplot(2,3,1)
                                maxD = max(disparities); minD = min(disparities);
                                plot(10.^disparities,marg_thr,'r'); hold on
                                text((10^maxD)/100,max(marg_thr)/2,sprintf('MAP: %2.f"',curr_est_max_thr))
                                text((10^maxD)/100,max(marg_thr)/3,['p(MAP+-',num2str(round(tolerance.*100)),'%) = ',sprintf('%2.1f%%',history(end,14))])
                                text((10^maxD)/100,max(marg_thr)/4,['CI_9_5% = [',sprintf('%2.f"',thr_CI(1)),'-',sprintf('%2.f"',thr_CI(2)),']'])
                                axis([10^minD 10^maxD 0 1.2*max(marg_thr)]);
                                plot([curr_est_max_thr curr_est_max_thr],[0 1],'r--')
                                xlabel('Thresholds (arcsec)')
                                title('Posterior for threshold, marginalized')
                                ylabel('Probability')
                                set(gca, 'XScale', 'log');
                                xticks([0.1 1 10 20 50 200 500 2000]);
                                xticklabels({'0.1' '1' '10' '20' '50' '200' '500' '2000'});
                                
                                subplot(2,3,2)
                                maxPS = max(psi.new_slopes); minPS = min(psi.new_slopes);
                                plot(psi.new_slopes,marg_pslo,'r'); hold on
                                text(maxPS/2,max(marg_pslo)/2,sprintf('MAP: %.2f',curr_est_max_pos_slo))
                                text(maxPS/2,max(marg_pslo)/3,['p(MAP+-',num2str(round(tolerance.*100)),'%) = ',sprintf('%2.1f%%',history(end,15))])
                                axis([minPS maxPS 0 1.2*max(marg_pslo)]);
                                plot([curr_est_max_pos_slo curr_est_max_pos_slo],[0 1],'r--')
                                xlabel('Positive slope')
                                title('Posterior for pos. slope, marginalized')
                                ylabel('Probability')
                                
                                subplot(2,3,3)
                                maxNS = max(psi.new_neg_slopes); minNS = min(psi.new_neg_slopes);
                                plot(psi.new_neg_slopes,marg_nslo,'r'); hold on
                                text(maxNS/2,max(marg_nslo)/2,sprintf('MAP: %.2f',curr_est_max_neg_slo))
                                text(maxNS/2,max(marg_nslo)/3,['p(MAP+-',num2str(round(tolerance.*100)),'%) = ',sprintf('%2.1f%%',history(end,16))])
                                axis([minNS maxNS 0 1.2*max(marg_nslo)]);
                                plot([curr_est_max_neg_slo curr_est_max_neg_slo],[0 1],'r--')
                                xlabel('Negative slope')
                                title('Posterior for neg. slope, marginalized')
                                ylabel('Probability')
                                                                
                                subplot(2,3,4);
                                for j=1:trial_nb
                                  if history(j,3)==1
                                    p1=plot(j,history(j,2),'og'); hold on
                                  else
                                    p2=plot(j,history(j,2),'xr'); hold on   
                                  end
                                end
                                p4 = plot(1:i,history(1:i,4),'r-'); 
                                xlabel('Trial')
                                ylabel('Disparity / threshold (arcsec)')
                                p3=plot(1:i,10.^history(1:i,9),'k-');
                                plot(1:i,10.^history(1:i,10),'k-');
                                p5=plot([1,i],[curr_est_max_thr curr_est_max_thr],'r--');
                                axis([1 trial_nb 0.1 10000])
                                set(gca, 'YScale', 'log')
                                yticks([0.1 1 10 20 50 200 500 2000]);
                                yticklabels({'0.1' '1' '10' '20' '50' '200' '500' '2000'});
                                title('Estimate by trial')
                                legend([p1, p2, p3, p4, p5],{'Correct','Incorrect','Grid range','MAP est.','Final est.'})
                                
                                subplot(2,3,5)
                                xxx = 10^minD:0.1:10^maxD;
                                plot(xxx, defineLikelihood_bell(g, curr_est_max_neg_slo, curr_est_max_pos_slo, delta, p, log10(xxx), log10(curr_est_max_thr), lapse),'r')
                                hold on
                                plot([curr_est_max_thr curr_est_max_thr],[0 1],'r--')
                                xlabel('Disparity (arcsec)')
                                ylabel('% CR')
                                text((10^maxD)/100, 0.35, sprintf('MAP: %d"',round(curr_est_max_thr)));
                                text((10^maxD)/100, 0.25,  sprintf('p(stereoblind) = %2.1f%%',100*sum(marg_thr((10.^disparities)>=1300))))
                                dispa = history(1:i,2); responses = history(1:i,3);
                                output = makeLevelEqualBoundsMean([dispa,responses],ceil(i/15));
                                scatter(output(:,1),output(:,2),output(:,3).*7,'ok')
                                set(gca, 'XScale', 'log');
                                xticks([0.1 1 10 20 50 200 500 2000]);
                                xticklabels({'0.1' '1' '10' '20' '50' '200' '500' '2000'});
                                title('Psychometric function')
                                
                                subplot(2,3,6);
                                plot(1:i,history(1:i,14),'r-'); hold on;
                                plot(1:i,history(1:i,15),'b-'); 
                                plot(1:i,history(1:i,16),'g-'); 
                                xlabel('Trial')
                                ylabel(['p(MAP+-',num2str(round(tolerance.*100)),'%) (%)'])
                                title('Confidence')
                                axis([1 trial_nb 0 100])
                                legend('threshold','pos. slope','neg. slope')
                       
                       end
                   end
     psi.threshold = curr_est_max_thr;
end