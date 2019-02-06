function [pp_mean,pc_mean] = ANALYSIS_1A_pcpp_mean(isubj,data)

	% This function calculates the mean for precahoust pump rt and prepump pump rt for each subject

	pcpp = zeros(length(data(isubj).PCPP),2)
	for i=1:length(data(isubj).PCPP)
		if data(isubj).PCPP(i) == 1
			pcpp(1,i) = data(isubj).rt(i)
			pcpp(2,i) = 0
		elseif data(isubj).PCPP(i) == 2
			pcpp(1,i) = 0
            		pcpp(2,i) = data(isubj).rt(i)
		end	
	end
	
	pp_mean = mean(nonzeros(pcpp(1,:)))
	pc_mean = mean(nonzeros(pcpp(2,:)))
