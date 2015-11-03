function inter = getIntervalCombinations(C,Er);
% find all possible intervals for max gap and pmax calculations
% note that we only want intervals with zero events for max gap
%
% 100414 pfs - modified old scripted code => function

		dis('finding all possible interval combinations...');
		clear inter logi
		n=0; % start with no events in each interval
		ii = 1; % starting point
		rr = 1; % running index of pairs
		while n<=(length(Er)-2) % all interval pairs
		%for n=1
			%try
			while ii<(length(Er)-n)
				inter.pair{rr} = [Er(ii) Er(ii+n+1)];
				inter.dEr{rr} = diff(inter.pair{rr});
				inter.n(rr) =  n;
				if 0 %%% figure out the time t of all events in the interval (including both endpoints)
					logi{rr} = Er>=inter.pair{rr}(1) & Er<=inter.pair{rr}(2);					
					all_t{rr} = data.days(logi{rr});
						[junk_r,tmp_c]=min(abs(all_t{rr} - 180));
						inter.min_t(rr) = all_t{rr}(tmp_c) / 365 + 0.415;
						inter.mean_t(rr) = mean(all_t{rr}) / 365 + 0.415;
					switch 1
					case 0 %  take the t closest to 180 (i.e. Dec 1).
						inter.t(rr) = inter.min_t(rr);
					case 1 %  take the t as average of all events in interval
						inter.t(rr) = inter.mean_t(rr);
					end						
				end
				ii=ii+1; % increment the loop counter
				rr=rr+1; % increment global counter
			end
			%catch
				%keyboard;
				ii=1; % re-set ii
				n=n+1; % increment n
			%end
		end
		inter.Er = Er; % record the observed event energies in this structure

		dis('done.');%return
