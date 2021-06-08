function psi = recomputeFromData4(psi, plotIt,titles)

%let's recompute everything from all data
psi.new_thresholds = unique([log10([1,10,100,1000,2000,10000]),linspace(log10(0.1),log10(3000),500),log10(psi.history(:,2))']); %range of possible values for thresholds T, in log10(arcsec)              
psi.new_slopes = linspace(0.2,2.4,30);
psi.new_neg_slopes = linspace(0,0.224,20);
psi.trial = size(psi.history,1); trial_nb=psi.trial;
storingDisp = nan(trial_nb,12); storingDisp(:,1:10)=psi.history;
%psi.new_thresholds = unique([psi1.new_thresholds,psi2.new_thresholds,psi.thresholdsIni]);
old_thresholds = psi.new_thresholds;
old_disparities = old_thresholds;
psi.new_disparities = [1 2]; %not important, will be simplified later
lapse = psi.lapse; g = psi.g; delta = psi.delta; p = psi.p;

               [psi.tt, psi.ss, psi.ll, psi.xx] = ndgrid(psi.new_thresholds, psi.new_slopes, psi.new_neg_slopes, psi.new_disparities);
                psi.prior = ones(size(psi.tt(:,:,:,1)))./numel(psi.tt(:,:,:,1)); % reinitialize a flat prior
                tt=psi.tt(:,:,:,1); ss = psi.ss(:,:,:,1); ll=psi.ll(:,:,:,1);
                ttt=tt(:);sss=ss(:);lll=ll(:);
                logUnit = 0.0212;                
                   for i=1:trial_nb % and update with the history of data
                      disp_j = psi.history(i,2); resp_j = psi.history(i,3);
                      psi.likelihoodCR = defineLikelihood_bell(psi.g, psi.ll(:,:,:,1), psi.ss(:,:,:,1), psi.delta, psi.p, log10(disp_j), psi.tt(:,:,:,1), psi.lapse); % psis for success
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
                      keyboard
                   if (plotIt==1) && (i==psi.trial)
                       % calculate stats for that trial (thresholds and marginal distributions)                  
                            weightedThres = (10.^tt).*psi.prior;    
                            sumWeightThres = sum(weightedThres(:));
                            psi.threshold = sumWeightThres;
                            curr_est_sum_thr = psi.threshold;
                            selected_thr=psi.prior(((psi.tt(:,1,1,1)>(log10(sumWeightThres)-logUnit))&(psi.tt(:,1,1,1)<(log10(sumWeightThres)+logUnit))),:,:,1);
                            if numel(selected_thr)==0; [~, selected_ind] = min(abs(psi.tt(:,1,1,1)-log10(sumWeightThres))); selected_thr = psi.prior(selected_ind,:,:,1); end
                            psi.pThreshold = sum(selected_thr(:)); % probability that this threshold is correct (cumulated prob sum of t+-0.0212 log unit ~+-5%)
                            
                            %marginalize distributions for stereoblindness calculation
                            marg_thr=squeeze(sum(sum(psi.prior(:,:,:,1),3),2));
                            storingDisp(i,11) = psi.pThreshold; 
                            storingDisp(i,7) = psi.threshold;
                            psi.stereoblind_prob = 100*sum(marg_thr((10.^psi.new_thresholds)>=psi.maxAllowerThreshold));  
                            
                             % max proba 
                            moving_aver = psi.prior+psi.prior([2:end,end],:,:)+psi.prior([1,1:(end-1)],:,:)+...
                                                   psi.prior(:,[2:end,end],:)+psi.prior(:,[1,1:(end-1)],:)+...
                                                   psi.prior(:,:,[2:end,end])+psi.prior(:,:,[1,1:(end-1)]);
                                               
                            [pCo2,idx2] = max(moving_aver(:)); 
                            storingDisp(i,12) = 100*pCo2/7; 
                            [pCo,idx3] = max(psi.prior(:)); 
                    end
                      
                      % if plotIt==1 && (mod(i,everyStep)==1 || i==psi.trial) 
                      if plotIt==1 && (i==psi.trial)
                           % max estimate
                             curr_est_max_thr = 10.^ttt(idx3);
                             curr_est_max_pos_slo = sss(idx3);
                             curr_est_max_neg_slo = lll(idx3);
                             
                           % estimate from max moving average (MA) proba  
                             curr_est_max_MA_thr = 10.^ttt(idx2);
                             curr_est_max_MA_pos_slo = sss(idx2);
                             curr_est_max_MA_neg_slo = lll(idx2);
                             
                            % estimate through weighted sum proba 
                            weightedSlope= psi.ss(:,:,:,1).*psi.prior;   
                            weightedNegSlo= psi.ll(:,:,:,1).*psi.prior;
                            curr_est_sum_pos_slo = sum(weightedSlope(:));
                            curr_est_sum_neg_slo = sum(weightedNegSlo(:));
                            

                           %     subplot(2,3,1)
%                                  hold off
%                                  plot(10.^old_thresholds(1:(end-1)),marg_thr(1:(end-1)),'r'); hold on
%                                  text(max(old_thresholds(1:(end-1)))/2,max(marg_thr(1:(end-1)))/2,sprintf('p(stereoblind) = %.2f%%',100*sum(marg_thr((10.^old_thresholds)>=1300))))
%                                 % plot([sim_threshold sim_threshold],[min(marg_thr), max(marg_thr(1:end))],'-k')
%                                 hold on
%                                 % plot([curr_est_sum_thr curr_est_sum_thr],[min(marg_thr), max(marg_thr(1:end))],'--r')
%                                 % plot(10.^[range(1) range(2)],[mean(marg_thr) mean(marg_thr)],'r-');
%                                  xlabel('thresholds (arcsec)')
%                                  title('Posterior for threshold, marginalized')
%                                  ylabel('Probability')
%                                  set(gca, 'XScale', 'log')
                                 
%                                 subplot(2,3,2)
%                                  hold off
%                                  plot(slopes,marg_slo,'r'); hold on
%                                  plot([sim_pos_slope sim_pos_slope],[min(marg_slo), max(marg_slo)],'-k')
%                                  plot([curr_est_sum_pos_slo curr_est_sum_pos_slo],[min(marg_slo), max(marg_slo)],'--r')
%                                  xlabel('Positive slopes')
%                                  title('Posterior for pos. slope, marginalized')
%                                  ylabel('Probability')
% 
%                                 subplot(2,3,3)
%                                 hold off
%                                  plot(neg_slopes,marg_lap,'r'); hold on
%                                  plot([sim_neg_slope sim_neg_slope],[min(marg_lap), max(marg_lap)],'-k')
%                                  plot([curr_est_sum_neg_slo curr_est_sum_neg_slo],[min(marg_lap), max(marg_lap)],'--r')
%                                  xlabel('Negative slopes')
%                                  title('Posterior for neg. slope, marginalized')
%                                  ylabel('Probability')
                             %   subplot(2,3,3)
%                                 axis([1 trial_nb 0 100])
%                                  hold on
%                                 plot(1:i,100*storingDisp(1:i,11),'r-');
%                                 plot(1:i,100*storingDisp(1:i,12),'b-');
%                                 xlabel('Trial')
%                                 ylabel('Certainty of estimate (%)')
%                                 legend('sum','max','Location','NorthWest')
                                
                             %   subplot(2,3,4)
                                 xxx = min(10.^old_disparities)/2:0.1:max(10.^old_disparities)*2;
                                % plot(xxx, defineLikelihood_bell(g, sim_neg_slope, sim_pos_slope, delta, p, log10(xxx), log10(sim_threshold), sim_lapse),'k')
                              %  hold off
                                  plot(xxx, defineLikelihood_bell(g, curr_est_max_neg_slo, curr_est_max_pos_slo, delta, p, log10(xxx), log10(curr_est_max_thr), lapse),'b')
                                  hold on
                                  plot(xxx, defineLikelihood_bell(g, curr_est_max_MA_neg_slo, curr_est_max_MA_pos_slo, delta, p, log10(xxx), log10(curr_est_max_MA_thr), lapse),'g')
                                  plot(xxx, defineLikelihood_bell(g, curr_est_sum_neg_slo, curr_est_sum_pos_slo, delta, p, log10(xxx), log10(curr_est_sum_thr), lapse),'r')
                                
                                 %plot([sim_threshold sim_threshold],[0 1],'k')
                                 plot([curr_est_max_thr curr_est_max_thr],[0 1],'b--')
                                 plot([curr_est_max_MA_thr curr_est_max_MA_thr],[0 1],'g--')
                                 plot([curr_est_sum_thr curr_est_sum_thr],[0 1],'r--')
                               %  legend('max', 'max MA', 'sum','Location','NorthWest')
                                 title(titles)
                                 xlabel('Log disparity (log arcsec)')
                                 ylabel('% CR')
                                %text(1, 0.7, sprintf('Trial %d',i))
                                 text(1, 0.65, sprintf('Max %d"',round(curr_est_max_thr)))
                                 text(1, 0.6, sprintf('Max MA %d"',round(curr_est_max_MA_thr)))
                                 text(1, 0.55, sprintf('Sum %d"',round(curr_est_sum_thr)))
                                 
                                 text(150, 0.65, sprintf('%1.2f%%',(pCo*100)))
                                 text(150, 0.6, sprintf('%1.2f%%',(pCo2*100/7)))
                                 text(150, 0.55, sprintf('%1.2f%%',(psi.pThreshold*100)))
                                 
                                 disparities = storingDisp(1:i,2); responses = storingDisp(1:i,3);
                                 output1 = makeLevelEqualBoundsMean([(disparities),responses],ceil(i/15));
                                 output2 = makeLevelEqualBoundsMean([(disparities),responses],ceil(i/15));
                                 output=[output1;output2]; %output(:,1)=10.^output(:,1);
                                 scatter(output(:,1),output(:,2),output(:,3).*7,'ok')
                                 axis([(10.^min(old_disparities)/2), (10.^max(old_disparities))*2, min([output(:,2)',g-0.1]), 1])
                                 set(gca, 'XScale', 'log')

%                               subplot(2,3,5);
%                               axis([1 trial_nb min(10.^storingDisp(1:i,11)) max(10.^storingDisp(1:i,12))])
%                               hold on
%                               plot(1:i,storingDisp(1:i,6),'r--');
%                               plot(1:i,storingDisp(1:i,3),'b--');     
%                               if storingDisp(i,2)==1
%                                 plot(i,storingDisp(i,1),'ok')
%                               else
%                                 plot(i,storingDisp(i,1),'xr')   
%                               end
%                               xlabel('Trial')
%                               ylabel('Disparity threshold')
%                              % plot([1 trial_nb], [sim_threshold sim_threshold],'k-')
%                                 plot(1:i,10.^storingDisp(1:i,11),'r-');
%                                 plot(1:i,10.^storingDisp(1:i,12),'r-');
%                                 legend('sum','max','correct/wrong','range');
%                                 set(gca, 'YScale', 'log')
% 
%                                subplot(2,3,6);
%                                axis([1 trial_nb min(storingDisp(1:i,13)) max(storingDisp(1:i,14))])
%                                hold on;
%                                plot(1:i,storingDisp(1:i,7),'r--');
%                              %  plot([1 trial_nb], [sim_pos_slope sim_pos_slope],'k-')
%                                plot(1:i,storingDisp(1:i,13),'r-');
%                                plot(1:i,storingDisp(1:i,14),'r-');
%                                xlabel('Trial')
%                                ylabel('Slope estimate')
%                                legend('sum','range');

                              %  drawnow

                       end
                      
                      
                      
                   end
            
                   
                   % calculate stats for that trial (thresholds and marginal distributions)                  
                            weightedThres = (10.^tt).*psi.prior;    
                            sumWeightThres = sum(weightedThres(:));
                            psi.threshold = sumWeightThres;
                            selected_thr=psi.prior(((psi.tt(:,1,1,1)>(log10(sumWeightThres)-logUnit))&(psi.tt(:,1,1,1)<(log10(sumWeightThres)+logUnit))),:,:,1);
                            if numel(selected_thr)==0; [~, selected_ind] = min(abs(psi.tt(:,1,1,1)-log10(sumWeightThres))); selected_thr = psi.prior(selected_ind,:,:,1); end
                            psi.pThreshold = sum(selected_thr(:)); % probability that this threshold is correct (cumulated prob sum of t+-0.0212 log unit ~+-5%)
                            
                            psi.threshold_max = curr_est_max_thr;
                            psi.threshold_max_MA = curr_est_max_MA_thr;
end