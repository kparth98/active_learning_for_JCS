classdef BanditArm
   properties
      pull_count
      mean_reward
      cumm_reward
      explore_bonus
      exp_bonus_const

      tx_beam
      rx_beam

      confidence
      sigma_h
   end

   methods
       function obj = BanditArm(tx_beam, rx_beam, confidence, sigma_h)
          obj.pull_count = 0;
          obj.mean_reward = 0;
          obj.cumm_reward = 0;
          obj.explore_bonus = Inf;
          N = length(tx_beam);
          obj.exp_bonus_const = log(pi*log2(N)/(3*confidence));
          
          if nargin>0
              obj.tx_beam = tx_beam;
              obj.rx_beam = rx_beam;
              obj.confidence = confidence;
              obj.sigma_h = sigma_h;
          else
              obj.tx_beam = 0;
              obj.rx_beam = 0;
              obj.confidence = 0;
              obj.sigma_h = 1;
          end
      end

      function ucb=UCB(obj)
        ucb = obj.mean_reward + obj.explore_bonus;
      end

      function lcb=LCB(obj)
          lcb = obj.mean_reward - obj.explore_bonus;
      end

      % function obj=pull_arm(obj,G,noise,t)
      %   reward = abs(obj.rx_beam'*(G*obj.tx_beam + noise));
      % 
      %   obj.cumm_reward = obj.cumm_reward + reward;
      %   obj.pull_count = obj.pull_count + 1;
      %   obj.mean_reward = obj.cumm_reward/obj.pull_count;
      % 
      %   % obj.mean_reward_squared = (reward^2 + obj.mean_reward_squared*obj.pull_count)/(obj.pull_count+1);
      %   % obj.var_reward = obj.mean_reward_squared - obj.mean_reward^2;
      % 
      %   obj.explore_bonus = obj.sigma_h*sqrt((4*log(pi*t^2/(3*obj.confidence)))/obj.pull_count);
      % end

      function obj=update_arm(obj,reward,t)
        obj.cumm_reward = obj.cumm_reward + reward;
        obj.pull_count = obj.pull_count + 1;
        obj.mean_reward = obj.cumm_reward/obj.pull_count;

        % obj.mean_reward_squared = (reward^2 + obj.mean_reward_squared*obj.pull_count)/(obj.pull_count+1);
        % obj.var_reward = obj.mean_reward_squared - obj.mean_reward^2;
      
        % obj.explore_bonus = obj.sigma_h*sqrt((4*log(pi*obj.pull_count^2/(3*obj.confidence)))/obj.pull_count);
        obj.explore_bonus = obj.sigma_h*sqrt(2*(obj.exp_bonus_const + 2*log(obj.pull_count))/obj.pull_count);
      end
   end
end