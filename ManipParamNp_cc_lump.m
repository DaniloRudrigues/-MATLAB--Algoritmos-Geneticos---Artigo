function [tags] = ManipParamNp_cc_lump(TagAnt,ac)
  
  npp = ac.npp;
  lump = ac.lump;
  npi = ac.npi;
  ncc = ac.ncc;
  Tc = ac.T;
  Mqp = ac.Qp;
  Mqi = ac.Qi;
  tags = [];
  x_min = 0.0001;
  
 opera = ac.opera;
 tempo = ac.time;

 if opera == 0            % ---> topado
%    for it = 1:ncc
%      prod_ac = 0;
%      inj_ac = 0;
%      text_p = 0;
%      text_i = 0;
%      p_f = 0;
%      if tempo == 1
%        if it == 1
%           time_days(1) = 0;
%        else
%           ivart = (npp + npi - 2)*ncc + it - 1;
%           time_days(it) = time_days(it-1) + Tc*365*TagAnt(ivart);
%        end 
%      end
%      if npp > 1            % ---> vazao produtores (absoluta)
%        p_i = (npp + npi - 2)*(it - 1) + 1;
%        p_f = p_i + npp - 2;
%        for ip = 1:(p_f - p_i +1)
% 	 if npi>1
%            cont_p = ip + (npp + npi)*(it - 1);
% 	 else
% 	   cont_p = ip + npp*(it - 1);
% 	 end
% 	 text_p = text_p + 1;
%          vazao_absoluta = TagAnt(p_i+ip-1)*Mqp;
%          tags(cont_p).name = sprintf('$GP%i_%i',it,text_p);
%          tags(cont_p).val = vazao_absoluta;
%          prod_ac = prod_ac + vazao_absoluta;
%          if tempo == 1
%            tags(cont_p).type = 1;                % ---> t_var
%            tags(cont_p).number = ip;             % (p_prod)
%            tags(cont_p).time = time_days(it);
%          end
%        end
%        vazao_ultimo_prod = max([Mqp*x_min (Mqp - prod_ac)]);
%        tags(cont_p+1).name = sprintf('$GP%i_%i',it,npp);
%        tags(cont_p+1).val = vazao_ultimo_prod;
%        if tempo == 1
%          tags(cont_p+1).type = 1;                 % ---> t_var
%          tags(cont_p+1).number = npp;             % (p_prod)
%          tags(cont_p+1).time = time_days(it);
%        end
%      end
%      
%      if npi > 1             % ---> vazao injetores (absoluta)
%        i_i = p_f + 1;
%        i_f = i_i + npi - 2;
%        for ii = 1:(i_f - i_i +1)
% 	 if npp>1
%            cont_i = ii + (npp + npi)*(it - 1) + npp;
% 	 else
% 	   cont_i = ii + npi*(it - 1);
% 	 end
%          text_i = text_i + 1;
%          vazao_absoluta = TagAnt(i_i+ii-1)*Mqi;
%          tags(cont_i).name = sprintf('$GI%i_%i',it,text_i);
%          tags(cont_i).val = vazao_absoluta;
%          inj_ac = inj_ac + vazao_absoluta;
%          if tempo == 1
%            tags(cont_i).type = 2;                  % ---> t_var
%            tags(cont_i).number = ii;               % (p_inj)
%            tags(cont_i).time = time_days(it);
%          end
%        end
%        vazao_ultimo_inj = max([Mqi*x_min (Mqi - inj_ac)]);
%        tags(cont_i+1).name = sprintf('$GI%i_%i',it,npi);
%        tags(cont_i+1).val = vazao_ultimo_inj;
%        if tempo == 1
%          tags(cont_i+1).type = 2;                   % ---> t_var
%          tags(cont_i+1).number = npi;               % (p_inj)
%          tags(cont_i+1).time = time_days(it);
%        end
%      end
%    end
 elseif opera == 1            % ---> n_topado
   for it = 1:ncc
     if tempo == 1
       if it == 1
          time_days(1) = 0;
       else
          ivart = (npp + npi)*ncc + it - 1;
          time_days(it) = time_days(it-1) + Tc*365*TagAnt(ivart);
       end
     end
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%% mexi aqui %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     for ip = 1:(lump*npp);     % ---> vazao produtores (absoluta)
       ivarp = ip + (lump*npp)*(it-1);
       vazao_absoluta = TagAnt(ivarp);
       tags(ivarp).name = sprintf('$GP%i_%i',it,ip);
       tags(ivarp).val = vazao_absoluta;
       if tempo == 1
         tags(ivarp).type = 1;      % ---> t_var
         tags(ivarp).number = ip;   % (p_prod)
         tags(ivarp).time = time_days(it);
       end
     end
%      for ip = 1:(npi);     % ---> vazao injetores (absoluta)
%        ivari = ip + (npp + npi)*(it-1) + (npp); 
%        vazao_absoluta = TagAnt(ivari)*Mqi;
%        tags(ivari).name = sprintf('$GI%i_%i',it,ip);
%        tags(ivari).val = vazao_absoluta;
%        if tempo == 1
%          tags(ivari).type = 2;      % ---> t_var
%          tags(ivari).number = ip;   % (p_inj)
%          tags(ivari).time = time_days(it);
%        end
     end
   end
 